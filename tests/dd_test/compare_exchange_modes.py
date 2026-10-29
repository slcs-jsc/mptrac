#!/usr/bin/env python3
"""Compare final MPTRAC binary snapshots by particle ID with bounded memory.

Each argument is a benchmark case directory (mode0 or mode1). The nine
column-major doubles are grouped into temporary ID buckets, so rank and
record order do not affect the comparison. Scratch space holds at most two
row-major copies of one snapshot and is removed on exit.
"""

import argparse
import array
import json
import math
import struct
import tempfile
from contextlib import ExitStack
from pathlib import Path


FIELDS = ("time", "pressure", "longitude", "latitude", "id",
          "tracer_0", "tracer_1", "tracer_2", "tracer_3")
NFIELD = len(FIELDS)
RECORD = struct.Struct("=9d")
CHUNK = 8192


def snapshot_files(case):
    rank_dirs = sorted((case / "data").glob("data.*"),
                       key=lambda p: int(p.name.split(".")[-1]))
    if not rank_dirs or [int(p.name.split(".")[-1]) for p in rank_dirs] != list(range(len(rank_dirs))):
        raise ValueError(f"Missing or nonconsecutive rank directories in {case}")
    names = [set(p.name for p in directory.glob("atm_*.bin"))
             for directory in rank_dirs]
    if any(group != names[0] for group in names) or not names[0]:
        raise ValueError(f"Snapshot names differ across ranks in {case}")
    stamp = max(names[0])
    return stamp, [directory / stamp for directory in rank_dirs]


def count_file(path):
    with path.open("rb") as stream:
        version, count = struct.unpack("=ii", stream.read(8))
        stream.seek(-4, 2)
        trailer, = struct.unpack("=i", stream.read(4))
    if version != 100 or count < 0 or trailer != 999 or \
            path.stat().st_size != 12 + RECORD.size * count:
        raise ValueError(f"Invalid MPTRAC binary snapshot: {path}")
    return count


def bucketize(paths, side, root, buckets, nonfinite):
    total = 0
    with ExitStack() as stack:
        outputs = [stack.enter_context((root / f"{side}-{i:03d}.bin").open("wb"))
                   for i in range(buckets)]
        for path in paths:
            count = count_file(path)
            total += count
            with path.open("rb") as stream:
                for start in range(0, count, CHUNK):
                    n = min(CHUNK, count - start)
                    columns = []
                    for field in range(NFIELD):
                        stream.seek(8 + 8 * (field * count + start))
                        values = array.array("d")
                        values.frombytes(stream.read(8 * n))
                        if len(values) != n:
                            raise ValueError(f"Truncated snapshot: {path}")
                        columns.append(values)
                    buffers = [bytearray() for _ in range(buckets)]
                    for row in zip(*columns):
                        pid = row[4]
                        if not math.isfinite(pid) or not pid.is_integer() or \
                                pid < 0 or pid > 2**53:
                            raise ValueError(f"Invalid particle ID {pid} in {path}")
                        for field, value in enumerate(row):
                            if math.isnan(value):
                                nonfinite["nan"][field] += 1
                            elif math.isinf(value):
                                nonfinite["inf"][field] += 1
                        buffers[int(pid) % buckets].extend(RECORD.pack(*row))
                    for output, buffer in zip(outputs, buffers):
                        if buffer:
                            output.write(buffer)
    return total


def compare_bucket(left_path, right_path, different_fields, maxima, max_bucket_bytes):
    if left_path.stat().st_size > max_bucket_bytes or \
            right_path.stat().st_size > max_bucket_bytes:
        raise ValueError("ID bucket exceeds memory cap; increase --buckets")
    left = left_path.read_bytes()
    lookup = {}
    duplicate_left = 0
    for offset in range(0, len(left), RECORD.size):
        pid, = struct.unpack_from("=d", left, offset + 4 * 8)
        key = int(pid)
        if key in lookup:
            duplicate_left += 1
        else:
            lookup[key] = offset
    duplicate_right = missing_left = different = 0
    bitwise = True
    seen_right = set()
    with right_path.open("rb") as stream:
        while data := stream.read(RECORD.size * CHUNK):
            for offset in range(0, len(data), RECORD.size):
                record = data[offset:offset + RECORD.size]
                pid, = struct.unpack_from("=d", record, 4 * 8)
                key = int(pid)
                if key in seen_right:
                    duplicate_right += 1
                    continue
                seen_right.add(key)
                old_offset = lookup.pop(key, None)
                if old_offset is None:
                    missing_left += 1
                    continue
                previous = left[old_offset:old_offset + RECORD.size]
                if previous != record:
                    bitwise = False
                old_values = RECORD.unpack(previous)
                new_values = RECORD.unpack(record)
                changed = False
                for field, (a, b) in enumerate(zip(old_values, new_values)):
                    if not math.isfinite(a) or not math.isfinite(b):
                        different_fields[field] += 1
                        changed = True
                    elif a != b:
                        different_fields[field] += 1
                        maxima[field] = max(maxima[field], abs(a - b))
                        changed = True
                different += changed
    return duplicate_left, duplicate_right, len(lookup), missing_left, different, bitwise


def compare(mode0, mode1, scratch, buckets, max_bucket_bytes):
    stamp0, files0 = snapshot_files(mode0)
    stamp1, files1 = snapshot_files(mode1)
    if stamp0 != stamp1:
        raise ValueError(f"Final snapshot timestamps differ: {stamp0}, {stamp1}")
    nonfinite = [{"nan": [0] * NFIELD, "inf": [0] * NFIELD}
                 for _ in range(2)]
    different_fields = [0] * NFIELD
    maxima = [0.0] * NFIELD
    totals = [0, 0, 0, 0, 0]
    bitwise = True
    with tempfile.TemporaryDirectory(prefix="mptrac-dd-compare-", dir=scratch) as tmp:
        root = Path(tmp)
        counts = [bucketize(files0, 0, root, buckets, nonfinite[0]),
                  bucketize(files1, 1, root, buckets, nonfinite[1])]
        for i in range(buckets):
            result = compare_bucket(root / f"0-{i:03d}.bin",
                                    root / f"1-{i:03d}.bin",
                                    different_fields, maxima, max_bucket_bytes)
            for j in range(5):
                totals[j] += result[j]
            bitwise &= result[5]
    duplicate0, duplicate1, only0, only1, different = totals
    same_ids = only0 == 0 and only1 == 0
    return {
        "snapshot": stamp0,
        "particle_count": {"mode0": counts[0], "mode1": counts[1]},
        "id_set_agreement": same_ids,
        "duplicate_ids": {"mode0": duplicate0, "mode1": duplicate1},
        "ids_only_in_mode0": only0,
        "ids_only_in_mode1": only1,
        "nonfinite": {f"mode{side}": {
            "nan": dict(zip(FIELDS, nonfinite[side]["nan"])),
            "inf": dict(zip(FIELDS, nonfinite[side]["inf"]))}
            for side in range(2)},
        "particles_with_differing_values": different,
        "differing_values_per_field": dict(zip(FIELDS, different_fields)),
        "max_absolute_difference_per_field": dict(zip(FIELDS, maxima)),
        "all_particle_records_bitwise_identical": same_ids and \
            duplicate0 == 0 and duplicate1 == 0 and bitwise,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode0", type=Path, help="mode0 case directory")
    parser.add_argument("mode1", type=Path, help="mode1 case directory")
    parser.add_argument("--scratch", type=Path, help="scratch filesystem with space for two final snapshots")
    parser.add_argument("--buckets", type=int, default=256,
                        help="temporary ID buckets, default 256 (roughly 32 MB per 100M bucket)")
    parser.add_argument("--max-bucket-mib", type=int, default=64,
                        help="maximum bytes per bucket loaded into memory, in MiB")
    args = parser.parse_args()
    if args.buckets < 1 or args.buckets > 512:
        parser.error("--buckets must be between 1 and 512")
    if args.max_bucket_mib < 1:
        parser.error("--max-bucket-mib must be positive")
    if args.scratch is not None and not args.scratch.is_dir():
        parser.error("--scratch must be an existing directory")
    result = compare(args.mode0, args.mode1, args.scratch, args.buckets,
                     args.max_bucket_mib * 2**20)
    print(json.dumps(result, indent=2, sort_keys=True))
    finite = all(not any(counts.values()) for mode in result["nonfinite"].values()
                 for counts in mode.values())
    unique = not any(result["duplicate_ids"].values())
    return 0 if result["id_set_agreement"] and unique and finite and \
        result["particles_with_differing_values"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
