#!/usr/bin/env python3
"""Aggregate exchange-mode campaign results.

Reads $DD_CAMPAIGN_ROOT/campaign_summary.csv (rows appended by
campaign_exchange_run.sbatch) and produces:

  campaign_aggregate.csv  mean and sample standard deviation of every timing
                          and count column, grouped by
                          (scenario, campaign_count, particles, dd_layout,
                          exchange_mode) across repeats.
  per_rank_timers.csv     raw per-rank TIMER_* measurements parsed from the
                          srun-labelled trac.log of every case directory
                          (scenario, campaign_count, total, layout, repeat,
                          mode, rank, timer, value).

Usage:
    ./campaign_exchange_report.py --root /path/to/campaign
"""
import argparse
import csv
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

KEYS = ("scenario", "campaign_count", "particles", "dd_layout", "exchange_mode")
TIMER_RE = re.compile(r"TIMER_([A-Z0-9_]+) = ([\d.eE+-]+) s")
RANK_RE = re.compile(r"^\s*(\d+):\s")


def numeric_columns(rows):
    columns = []
    for name in rows[0]:
        if name in KEYS or name in ("repeat", "run_id"):
            continue
        values = [row[name] for row in rows if row[name] not in ("", None)]
        if not values:
            continue
        try:
            [float(value) for value in values]
        except ValueError:
            continue
        columns.append(name)
    return columns


def aggregate(root):
    summary = root / "campaign_summary.csv"
    if not summary.is_file():
        raise SystemExit(f"[ERROR] Missing {summary}; run campaign cells first")
    with open(summary, newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise SystemExit(f"[ERROR] Empty {summary}")
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[key] for key in KEYS)].append(row)
    measure = numeric_columns(rows)
    out = root / "campaign_aggregate.csv"
    with open(out, "w", newline="") as stream:
        fieldnames = list(KEYS) + ["n_runs"]
        for name in measure:
            fieldnames += [f"{name}_mean", f"{name}_stdev"]
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for key in sorted(groups, key=lambda k: (k[0], int(k[1]), int(k[2]),
                                                 k[3], int(k[4]))):
            members = groups[key]
            record = dict(zip(KEYS, key))
            record["n_runs"] = len(members)
            for name in measure:
                values = [float(row[name]) for row in members
                          if row[name] not in ("", None)]
                record[f"{name}_mean"] = (f"{statistics.fmean(values):.6g}"
                                          if values else "")
                record[f"{name}_stdev"] = (f"{statistics.stdev(values):.6g}"
                                           if len(values) > 1 else "")
            writer.writerow(record)
    print(f"[INFO] Wrote {out} ({len(groups)} groups)")


def per_rank_timers(root):
    out = root / "per_rank_timers.csv"
    fields = ["scenario", "campaign_count", "total", "layout", "repeat",
              "mode", "rank", "timer", "seconds"]
    count = 0
    with open(out, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        for scen_dir in sorted((root / "runs").glob("*")):
            if not scen_dir.is_dir():
                continue
            for cell_dir in sorted(scen_dir.glob("c*")):
                count_name = cell_dir.name[1:]
                for total_dir in sorted(cell_dir.glob("n*")):
                    for layout_dir in sorted(total_dir.glob("layout*")):
                        for repeat_dir in sorted(layout_dir.glob("repeat*")):
                            for mode_dir in sorted(repeat_dir.glob("mode*")):
                                log = mode_dir / "trac.log"
                                if not log.is_file():
                                    continue
                                for line in log.read_text(
                                        errors="replace").splitlines():
                                    timer = TIMER_RE.search(line)
                                    if not timer:
                                        continue
                                    rank = RANK_RE.match(line)
                                    writer.writerow([
                                        scen_dir.name, count_name,
                                        total_dir.name[1:],
                                        layout_dir.name[len("layout"):],
                                        repeat_dir.name[len("repeat"):],
                                        mode_dir.name[len("mode"):],
                                        rank.group(1) if rank else "na",
                                        timer.group(1), timer.group(2)])
                                    count += 1
    print(f"[INFO] Wrote {out} ({count} per-rank timer samples)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path,
                        help="Campaign root directory")
    args = parser.parse_args()
    root = args.root.resolve()
    aggregate(root)
    per_rank_timers(root)


if __name__ == "__main__":
    sys.exit(main())
