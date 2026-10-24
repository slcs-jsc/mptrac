#!/usr/bin/env python3
"""Validate binary DD snapshots and summarize existing MPTRAC timers."""

import argparse
import array
import csv
import math
import re
import struct
import sys
from pathlib import Path


FIELDS = (
    "particles", "mpi_ranks", "gpus", "dd_layout", "exchange_mode",
    "sort_dt", "simulation_hours", "total_runtime", "launcher_wall_time",
    "dd_runtime",
    "classification_time", "sorting_time", "compaction_time", "packing_time",
    "mpi_time", "gpu_transfer_time", "gpu_d2h_time",
    "pack_transfer_phase_time", "unpack_h2d_phase_time",
    "migrated_particles", "migration_fraction", "peak_particles_per_rank",
    "initial_particles_per_rank", "estimated_particle_memory_gib_per_rank",
    "sent_particle_events", "repeat",
)


def header(path, nq):
    with path.open("rb") as stream:
        version, count = struct.unpack("=ii", stream.read(8))
    if version != 100 or count < 0:
        raise ValueError(f"Invalid MPTRAC binary header: {path}")
    expected = 12 + 8 * (4 + nq) * count
    if path.stat().st_size != expected:
        raise ValueError(f"Unexpected binary file size: {path}")
    return count


def particle_ids(path, count):
    # The first quantity is idx. Skip the four geolocation arrays and read
    # only idx in bounded chunks; no complete particle output enters memory.
    with path.open("rb") as stream:
        stream.seek(8 + 32 * count)
        remaining = count
        while remaining:
            nread = min(remaining, 131072)
            values = array.array("d")
            values.frombytes(stream.read(8 * nread))
            if len(values) != nread:
                raise ValueError(f"Truncated particle IDs: {path}")
            yield from values
            remaining -= nread


def timers(path):
    result = {}
    peaks = []
    sent = 0
    for line in path.read_text(errors="replace").splitlines():
        match = re.search(r"TIMER_([A-Z0-9_]+) = ([\d.eE+-]+) s", line)
        if match:
            name, value = match.group(1), float(match.group(2))
            detail = re.search(r"mean=\s*([\d.eE+-]+) s, max=.*?n=\s*(\d+)\)", line)
            if detail:
                value = float(detail.group(1)) * int(detail.group(2))
            result[name] = max(result.get(name, 0.0), value)
        match = re.search(r"SIZE_NP_PEAK = (\d+)", line)
        if match:
            peaks.append(int(match.group(1)))
        match = re.search(r"DD_SENT_PARTICLES = (\d+)", line)
        if match:
            sent += int(match.group(1))
    return result, max(peaks) if peaks else None, sent


def validate_snapshots(root, ranks, per_rank, stride, nq):
    expected = ranks * per_rank
    rank_dirs = [root / "data" / f"data.{rank}" for rank in range(ranks)]
    if any(not directory.is_dir() for directory in rank_dirs):
        raise ValueError("Missing rank directory")
    for directory in rank_dirs:
        if header(directory / "init.bin", nq) != per_rank:
            raise ValueError(f"Wrong initial particle count in {directory}")
        if list(directory.glob("atm_*.tab")):
            raise ValueError(f"ASCII particle output in benchmark: {directory}")

    stamps = sorted({path.name for directory in rank_dirs
                     for path in directory.glob("atm_*.bin")})
    if len(stamps) < 2:
        raise ValueError("Expected initial and final binary snapshots")
    migrated = 0
    for stamp in stamps:
        seen = bytearray((expected + 7) // 8)
        count_all = 0
        for rank, directory in enumerate(rank_dirs):
            path = directory / stamp
            count = header(path, nq)
            count_all += count
            for value in particle_ids(path, count):
                if not math.isfinite(value) or int(value) != value:
                    raise ValueError(f"Invalid particle ID in {path}")
                pid = int(value)
                origin, local = divmod(pid, stride)
                if origin >= ranks or local >= per_rank or origin < 0:
                    raise ValueError(f"Out-of-range particle ID {pid} in {path}")
                ordinal = origin * per_rank + local
                byte, bit = divmod(ordinal, 8)
                mask = 1 << bit
                if seen[byte] & mask:
                    raise ValueError(f"Duplicate particle ID {pid} at {stamp}")
                seen[byte] |= mask
                if stamp == stamps[-1] and origin != rank:
                    migrated += 1
        if count_all != expected:
            raise ValueError(f"Lost particles at {stamp}: {count_all} != {expected}")
    return migrated, len(stamps)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--ranks", type=int, required=True)
    parser.add_argument("--gpus", type=int, required=True)
    parser.add_argument("--lon-domains", type=int, required=True)
    parser.add_argument("--lat-domains", type=int, required=True)
    parser.add_argument("--per-rank", type=int, required=True)
    parser.add_argument("--id-stride", type=int, required=True)
    parser.add_argument("--capacity", type=int, required=True)
    parser.add_argument("--nq", type=int, required=True)
    parser.add_argument("--mode", type=int, required=True)
    parser.add_argument("--sort-dt", type=int, required=True)
    parser.add_argument("--hours", type=int, required=True)
    parser.add_argument("--repeat", type=int, required=True)
    parser.add_argument("--header", action="store_true")
    args = parser.parse_args()

    migrated, snapshots = validate_snapshots(
        args.case, args.ranks, args.per_rank, args.id_stride, args.nq)
    if migrated == 0:
        raise ValueError("No particles migrated between initial and final ranks")
    measured, peak, sent = timers(args.log)
    if "GROUP_DD" not in measured or "TOTAL" not in measured:
        raise ValueError("Missing DD or total MPTRAC timer output")
    if sent < migrated:
        raise ValueError("Missing or inconsistent DD send-count diagnostics")
    total = args.per_rank * args.ranks
    wall_file = args.case / "launcher_wall_seconds.txt"
    # atm_t, cache_t, and the three DD helper arrays. Excludes meteorology,
    # MPI buffers, allocator overhead, and possible host/device duplicates.
    bytes_per_slot = 8 * (4 + args.nq) + 60 + 20
    value = lambda name: measured.get(name, "")
    row = {
        "particles": total,
        "mpi_ranks": args.ranks,
        "gpus": args.gpus,
        "dd_layout": f"{args.lon_domains}x{args.lat_domains}",
        "exchange_mode": args.mode,
        "sort_dt": args.sort_dt,
        "simulation_hours": args.hours,
        "total_runtime": value("TOTAL"),
        "launcher_wall_time": float(wall_file.read_text()) if wall_file.exists() else "",
        "dd_runtime": value("GROUP_DD"),
        "classification_time": value("DD_GPU_CLASSIFY") or value("DD_COUNT_OUTGOING"),
        "sorting_time": value("DD_SORT") if args.mode == 0 else "",
        "compaction_time": value("DD_GPU_COMPACT") if args.mode == 1 else "",
        "packing_time": value("DD_GPU_PACK") if args.mode == 1 else "",
        # This timer includes host-side preparation in addition to MPI calls.
        "mpi_time": value("DD_COMMUNICATE_PARTICLES"),
        # A complete GPU transfer total is not timed independently.
        "gpu_transfer_time": "",
        "gpu_d2h_time": value("DD_GPU_SEND_D2H"),
        "pack_transfer_phase_time": value("DD_ATM2PARTICLES"),
        "unpack_h2d_phase_time": value("DD_PARTICLES2ATM"),
        "migrated_particles": migrated,
        "migration_fraction": migrated / total,
        "peak_particles_per_rank": peak if peak is not None else "",
        "initial_particles_per_rank": args.per_rank,
        "estimated_particle_memory_gib_per_rank":
            round(bytes_per_slot * args.capacity / 2**30, 6),
        "sent_particle_events": sent,
        "repeat": args.repeat,
    }
    for key in ("total_runtime", "launcher_wall_time", "dd_runtime", "classification_time",
                "sorting_time", "compaction_time", "packing_time",
                "mpi_time", "gpu_d2h_time", "pack_transfer_phase_time",
                "unpack_h2d_phase_time"):
        if row[key] != "":
            row[key] = f"{row[key]:.9g}"
    row["migration_fraction"] = f"{row['migration_fraction']:.9g}"
    writer = csv.DictWriter(sys.stdout, fieldnames=FIELDS, lineterminator="\n")
    if args.header:
        writer.writeheader()
    writer.writerow(row)
    print(f"Validated {total} IDs in {snapshots} binary snapshots; "
          f"{migrated} changed owner, {sent} send events", file=sys.stderr)


if __name__ == "__main__":
    main()
