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
    "particles", "nodes", "mpi_ranks", "gpus", "ranks_per_node",
    "dd_layout", "exchange_mode",
    "sort_dt", "simulation_hours", "total_runtime", "launcher_wall_time",
    "dd_runtime",
    "classification_time", "sorting_time", "compaction_time", "packing_time",
    "assign_time", "push_time",
    "mpi_time", "gpu_transfer_time", "gpu_d2h_time",
    "pack_transfer_phase_time", "unpack_h2d_phase_time",
    "migrated_particles", "migration_fraction", "peak_particles_per_rank",
    "initial_particles_per_rank", "estimated_particle_memory_gib_per_rank",
    "sent_particle_events", "intra_node_particle_events",
    "inter_node_particle_events", "intra_node_bytes", "inter_node_bytes",
    "inter_node_event_fraction", "repeat",
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


def check_finite_coordinates(path, count):
    with path.open("rb") as stream:
        stream.seek(8)
        for field in ("time", "pressure", "longitude", "latitude"):
            remaining = count
            while remaining:
                nread = min(remaining, 131072)
                values = array.array("d")
                values.frombytes(stream.read(8 * nread))
                if len(values) != nread or any(not math.isfinite(value)
                                                for value in values):
                    raise ValueError(f"Non-finite {field} in {path}")
                remaining -= nread


def timers(path):
    result = {}
    peaks = []
    sent = 0
    mapping = {}
    communication = {}
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
        match = re.search(
            r"(?:DD_SENT_PARTICLES =\s*|DD_SEND rank=\d+ particles=)(\d+)",
            line)
        if match:
            sent += int(match.group(1))
        match = re.search(
            r"(?:DD_BENCH_NODE|DD_NODE) rank=(\d+) leader=(\d+) local_rank=(\d+) "
            r"host=(\S+) gpu=(\S+)", line)
        if match:
            rank, leader, local_rank = map(int, match.group(1, 2, 3))
            if rank in mapping:
                raise ValueError(f"Duplicate node mapping for rank {rank}")
            mapping[rank] = (leader, local_rank, match.group(4), match.group(5))
        match = re.search(
            r"(?:DD_BENCH_COMM|DD_COMM) rank=(\d+) sent=(\d+) intra=(\d+) "
            r"inter=(\d+) intra_bytes=(\d+) inter_bytes=(\d+) "
            r"particle_bytes=(\d+)", line)
        if match:
            values = tuple(map(int, match.groups()))
            if values[0] in communication:
                raise ValueError(f"Duplicate communication summary for rank {values[0]}")
            communication[values[0]] = values[1:]
    return result, max(peaks) if peaks else None, sent, mapping, communication


def validate_communication(mapping, communication, ranks, nodes, ranks_per_node,
                           gpus, case):
    if set(mapping) != set(range(ranks)) or set(communication) != set(range(ranks)):
        raise ValueError("Missing DD rank placement or communication diagnostics; "
                         "compile MPTRAC with LOGLEV >= 2")
    by_node = {}
    for rank, (leader, local_rank, host, gpu) in mapping.items():
        by_node.setdefault(leader, []).append((rank, local_rank, host, gpu))
    if len(by_node) != nodes or any(len(group) != ranks_per_node
                                   for group in by_node.values()):
        raise ValueError(f"Actual MPI placement differs from {nodes} nodes × "
                         f"{ranks_per_node} ranks/node")
    for leader, group in by_node.items():
        if leader not in {item[0] for item in group} or \
                set(item[1] for item in group) != set(range(ranks_per_node)) or \
                len({item[2] for item in group}) != 1:
            raise ValueError(f"Inconsistent shared-memory mapping for node {leader}")
        visible = [item[3] for item in group]
        if gpus == ranks and all(gpu != "-" for gpu in visible):
            if any(len(gpu.split(",")) != 1 for gpu in visible) or \
                    len(set(visible)) != len(visible):
                raise ValueError(f"GPU assignment overlaps on node {leader}")
    with (case / "rank_mapping.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("rank", "node_leader_rank", "local_rank", "host",
                         "cuda_visible_devices"))
        for rank in sorted(mapping):
            writer.writerow((rank, *mapping[rank]))
    sums = [sum(communication[rank][i] for rank in range(ranks))
            for i in range(5)]
    sent, intra, inter, intra_bytes, inter_bytes = sums
    particle_sizes = {communication[rank][5] for rank in range(ranks)}
    if len(particle_sizes) != 1 or sent != intra + inter or \
            intra_bytes != intra * next(iter(particle_sizes)) or \
            inter_bytes != inter * next(iter(particle_sizes)):
        raise ValueError("Inconsistent MPI particle counts or payload bytes")
    if nodes == 1 and inter != 0:
        raise ValueError("Single-node run reported inter-node particle transfers")
    return sent, intra, inter, intra_bytes, inter_bytes


def initial_counts(init_root, ranks, nq):
    return [header(Path(init_root) / f"data.{rank}" / "init.bin", nq)
            for rank in range(ranks)]


def build_origin_map(init_root, ranks, nq):
    counts = initial_counts(init_root, ranks, nq)
    total = sum(counts)
    origins = array.array("h", [-1]) * total
    for rank in range(ranks):
        path = Path(init_root) / f"data.{rank}" / "init.bin"
        for value in particle_ids(path, counts[rank]):
            pid = int(value)
            if not math.isfinite(value) or int(value) != value or \
                    pid < 0 or pid >= total or origins[pid] != -1:
                raise ValueError(f"Invalid or duplicate initial particle ID "
                                 f"{value} in {path}")
            origins[pid] = rank
    return origins, counts


def validate_snapshots(root, ranks, per_rank, stride, nq, finite_coordinates=False,
                       init_source=None):
    origins = None
    expected = ranks * per_rank
    if init_source is not None:
        origins, counts = build_origin_map(init_source, ranks, nq)
        expected = sum(counts)
    rank_dirs = [root / "data" / f"data.{rank}" for rank in range(ranks)]
    if any(not directory.is_dir() for directory in rank_dirs):
        raise ValueError("Missing rank directory")
    for rank, directory in enumerate(rank_dirs):
        count = header(directory / "init.bin", nq)
        if origins is None:
            if count != per_rank:
                raise ValueError(f"Wrong initial particle count in {directory}")
        elif count != counts[rank]:
            raise ValueError(f"Initial particle count differs from the init "
                             f"source in {directory}")
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
            if finite_coordinates:
                check_finite_coordinates(path, count)
            count_all += count
            for value in particle_ids(path, count):
                if not math.isfinite(value) or int(value) != value:
                    raise ValueError(f"Invalid particle ID in {path}")
                pid = int(value)
                if origins is None:
                    origin, local = divmod(pid, stride)
                    if origin >= ranks or local >= per_rank or origin < 0:
                        raise ValueError(f"Out-of-range particle ID {pid} in {path}")
                    ordinal = origin * per_rank + local
                else:
                    origin = origins[pid] if 0 <= pid < expected else -1
                    if origin < 0:
                        raise ValueError(f"Out-of-range particle ID {pid} in {path}")
                    ordinal = pid
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
    parser.add_argument("--nodes", type=int, required=True)
    parser.add_argument("--ranks-per-node", type=int, required=True)
    parser.add_argument("--gpus", type=int, required=True)
    parser.add_argument("--lon-domains", type=int, required=True)
    parser.add_argument("--lat-domains", type=int, required=True)
    parser.add_argument("--per-rank", type=int, required=True)
    parser.add_argument("--id-stride", type=int, required=True)
    parser.add_argument("--init-source", type=Path)
    parser.add_argument("--capacity", type=int, required=True)
    parser.add_argument("--nq", type=int, required=True)
    parser.add_argument("--mode", type=int, required=True)
    parser.add_argument("--sort-dt", type=int, required=True)
    parser.add_argument("--hours", type=int, required=True)
    parser.add_argument("--repeat", type=int, required=True)
    parser.add_argument("--allow-zero-migration", action="store_true")
    parser.add_argument("--check-finite-coordinates", action="store_true")
    parser.add_argument("--header", action="store_true")
    args = parser.parse_args()

    migrated, snapshots = validate_snapshots(
        args.case, args.ranks, args.per_rank, args.id_stride, args.nq,
        args.check_finite_coordinates, args.init_source)
    if migrated == 0 and not args.allow_zero_migration:
        raise ValueError("No particles migrated between initial and final ranks")
    measured, peak, old_sent, mapping, communication = timers(args.log)
    sent, intra, inter, intra_bytes, inter_bytes = validate_communication(
        mapping, communication, args.ranks, args.nodes, args.ranks_per_node,
        args.gpus, args.case)
    if "GROUP_DD" not in measured or "TOTAL" not in measured:
        raise ValueError("Missing DD or total MPTRAC timer output")
    if sent < migrated or (old_sent and old_sent != sent):
        raise ValueError("Missing or inconsistent DD send-count diagnostics")
    if args.init_source is not None:
        total = sum(initial_counts(args.init_source, args.ranks, args.nq))
    else:
        total = args.per_rank * args.ranks
    wall_file = args.case / "launcher_wall_seconds.txt"
    # atm_t, cache_t, and the three DD helper arrays. Excludes meteorology,
    # MPI buffers, allocator overhead, and possible host/device duplicates.
    bytes_per_slot = 8 * (4 + args.nq) + 60 + 20
    value = lambda name: measured.get(name, "")
    row = {
        "particles": total,
        "nodes": args.nodes,
        "mpi_ranks": args.ranks,
        "gpus": args.gpus,
        "ranks_per_node": args.ranks_per_node,
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
        "assign_time": value("DD_ASSIGN_SUBDOMAINS"),
        "push_time": value("DD_PUSH"),
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
        "intra_node_particle_events": intra,
        "inter_node_particle_events": inter,
        "intra_node_bytes": intra_bytes,
        "inter_node_bytes": inter_bytes,
        "inter_node_event_fraction": inter / sent if sent else 0,
        "repeat": args.repeat,
    }
    for key in ("total_runtime", "launcher_wall_time", "dd_runtime", "classification_time",
                "sorting_time", "compaction_time", "packing_time", "assign_time",
                "push_time",
                "mpi_time", "gpu_d2h_time", "pack_transfer_phase_time",
                "unpack_h2d_phase_time"):
        if row[key] != "":
            row[key] = f"{row[key]:.9g}"
    row["migration_fraction"] = f"{row['migration_fraction']:.9g}"
    row["inter_node_event_fraction"] = f"{row['inter_node_event_fraction']:.9g}"
    writer = csv.DictWriter(sys.stdout, fieldnames=FIELDS, lineterminator="\n")
    if args.header:
        writer.writeheader()
    writer.writerow(row)
    print(f"Validated {total} IDs in {snapshots} binary snapshots; "
          f"{migrated} changed owner, {sent} send events "
          f"({intra} intra-node, {inter} inter-node)", file=sys.stderr)


if __name__ == "__main__":
    main()
