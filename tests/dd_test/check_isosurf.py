#!/usr/bin/env python3
"""Compare isosurface DD trajectories with independent runs without DD."""

import math
import sys
from pathlib import Path


def snapshots(root):
    result = {}
    for rank_dir in sorted(root.glob("data.*")):
        rank = int(rank_dir.name.split(".")[-1])
        for path in rank_dir.glob("atm_*.tab"):
            stamp = path.stem
            particles = result.setdefault(stamp, {})
            for line in path.read_text().splitlines():
                if not line or line.startswith("#"):
                    continue
                fields = line.split()
                pid = int(float(fields[4]))
                if pid in particles:
                    raise AssertionError(f"duplicate particle {pid} in {stamp}")
                particles[pid] = (rank, [float(value) for value in fields])
    return result


def compare(dd_root, serial_root, label):
    dd = snapshots(dd_root)
    serial = snapshots(serial_root)
    assert dd.keys() == serial.keys(), f"{label}: output times differ"
    assert len(dd) >= 2, f"{label}: need multiple output times"
    stamps = sorted(dd)
    initial_rank = {pid: rank for pid, (rank, _) in dd[stamps[0]].items()}
    assert len(initial_rank) == 48, f"{label}: expected 48 initial particles"
    iso_values = {row[9] for _, row in dd[stamps[0]].values()}
    assert len(iso_values) > 1, f"{label}: reference values do not vary"

    for stamp in stamps:
        assert dd[stamp].keys() == serial[stamp].keys() == initial_rank.keys(), (
            f"{label}: lost or extra particles at {stamp}"
        )
        for pid, (_, row) in dd[stamp].items():
            reference = serial[stamp][pid][1]
            # Columns: time, altitude, longitude, latitude, idx, zeta, mass,
            # current_subdomain, target_subdomain, iso_var.
            for col in (0, 1, 2, 3, 5, 6, 9):
                if not math.isclose(row[col], reference[col], rel_tol=1e-9,
                                    abs_tol=1e-6):
                    raise AssertionError(
                        f"{label}: particle {pid}, {stamp}, column {col}: "
                        f"DD={row[col]}, serial={reference[col]}"
                    )
            assert row[9] == dd[stamps[0]][pid][1][9], (
                f"{label}: iso_var changed for particle {pid}"
            )

    migrated = sum(initial_rank[pid] != rank
                   for pid, (rank, _) in dd[stamps[-1]].items())
    assert migrated > 0, f"{label}: no particles migrated between MPI ranks"
    print(f"{label}: {len(initial_rank)} particles, {len(stamps)} snapshots, "
          f"{migrated} migrated; iso_var and trajectories match without DD")


if __name__ == "__main__":
    compare(Path(sys.argv[1]), Path(sys.argv[3]), "DD_EXCHANGE_MODE=0")
    compare(Path(sys.argv[2]), Path(sys.argv[3]), "DD_EXCHANGE_MODE=1")
