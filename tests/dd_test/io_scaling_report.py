#!/usr/bin/env python3
"""Summarize the 137-level DD NetCDF I/O scaling runs."""

import argparse
import csv
import html
import re
import statistics
from pathlib import Path


LAYOUTS = ("1x4", "2x2", "4x1")
FIELDS = (
    "resolution", "nx", "ny", "nz", "layout", "repeat", "status",
    "io_trace_status", "gpu_sample_status",
    "particles", "total_seconds", "input_seconds", "dd_seconds",
    "core_seconds_max_rank", "halo_seconds_max_rank", "close_seconds_max_rank",
    "requested_bytes_all_ranks", "requested_gib_per_input_second",
    "core_calls_all_ranks", "halo_calls_all_ranks", "core_shape", "halo_shape",
    "peak_host_gib_per_rank", "peak_gpu_gib", "wind_files_gib",
    "sent_particle_events", "migration_fraction",
)


def input_timer(path):
    if not path.exists():
        return ""
    values = [float(x) for x in re.findall(
        r"TIMER_GROUP_INPUT =\s*([\d.eE+-]+) s", path.read_text(errors="replace"))]
    return max(values) if values else ""


def rank_timers(path):
    timers = {}
    if not path.exists():
        return timers
    for line in path.read_text(errors="replace").splitlines():
        match = re.match(r"\s*(\d+): TIMER_(GROUP_INPUT|GROUP_DD|TOTAL) =\s*([\d.eE+-]+) s", line)
        if match:
            rank, name, value = match.groups()
            timers.setdefault(int(rank), {})[name.lower() + "_seconds"] = float(value)
    return timers


def traces(path):
    files = sorted(path.glob("rank_*.csv"))
    if len(files) != 4:
        return {}
    per_rank = []
    for file in files:
        with file.open(newline="") as stream:
            records = list(csv.DictReader(stream))
        if len(records) != 7:
            return {}
        per_rank.append(records)
    sums = lambda key: [sum(float(record[key]) for record in records)
                        for records in per_rank]
    halo_bytes = int(sum(sums("halo_bytes")))
    shapes = lambda prefix: ";".join(sorted({"x".join(
        record[f"{prefix}_{dim}"] for dim in ("lev", "lat", "lon"))
        for records in per_rank for record in records
        if record[f"{prefix}_lev"] != "0"}))
    return {
        "core_seconds_max_rank": max(sums("core_seconds")),
        "halo_seconds_max_rank": max(sums("halo_seconds")) if halo_bytes else 0,
        "close_seconds_max_rank": max(sums("close_seconds")),
        "requested_bytes_all_ranks": int(sum(sums("core_bytes"))) + halo_bytes,
        "core_calls_all_ranks": int(sum(sums("core_calls"))),
        "halo_calls_all_ranks": int(sum(sums("halo_calls"))) if halo_bytes else 0,
        "core_shape": shapes("core"),
        "halo_shape": shapes("halo"),
        "peak_host_gib_per_rank": max(float(record["peak_rss_kib"])
            for records in per_rank for record in records) / 2**20,
    }


def rank_traces(path, resolution, layout, repeat):
    rows = []
    for file in sorted(path.glob("rank_*.csv")):
        with file.open(newline="") as stream:
            records = list(csv.DictReader(stream))
        if len(records) != 7:
            continue
        row = {"resolution": resolution, "layout": layout, "repeat": repeat,
               "rank": int(file.stem.split("_")[-1])}
        for field in ("core_seconds", "halo_seconds", "close_seconds",
                      "core_calls", "halo_calls", "core_bytes", "halo_bytes"):
            row[field] = sum(float(record[field]) for record in records)
        row["peak_host_gib"] = max(float(record["peak_rss_kib"])
                                    for record in records) / 2**20
        rows.append(row)
    return rows


def gpu_peak(path):
    if not path.exists():
        return ""
    with path.open(newline="") as stream:
        values = [float(row["memory_mib"]) for row in csv.DictReader(stream)
                  if row.get("memory_mib", "").strip().isdigit()]
    return max(values) / 1024 if values else ""


def plot(rows, field, ylabel, path):
    series = []
    for layout in LAYOUTS:
        points = []
        for resolution in (1.0, 0.5, 0.25, 0.125):
            values = [float(row[field]) for row in rows
                      if row["status"] == "complete" and
                      row["layout"] == layout and
                      float(row["resolution"]) == resolution and
                      row[field] != ""]
            if values:
                points.append((1 / resolution, statistics.median(values)))
        if points:
            series.append((layout, points))
    if not series:
        return ""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        svg = path.with_suffix(".svg")
        maximum = max(value for _, points in series for _, value in points) * 1.1
        maximum = maximum or 1
        colors = ("#1f77b4", "#d95f02", "#2ca02c")
        lines = ['<svg xmlns="http://www.w3.org/2000/svg" width="720" height="430" viewBox="0 0 720 430">',
                 '<rect width="720" height="430" fill="white"/>',
                 '<path d="M75 30V360H680" fill="none" stroke="black"/>',
                 f'<text x="20" y="20" font-size="16">{html.escape(ylabel)}</text>']
        for tick, label in zip((1, 2, 4, 8), ("1°", "0.5°", "0.25°", "0.125°")):
            x = 75 + (tick - 1) * 85
            lines.append(f'<text x="{x-16}" y="382" font-size="13">{label}</text>')
        for value in (0, maximum / 2, maximum):
            y = 360 - 320 * value / maximum
            lines.append(f'<text x="8" y="{y+5:.1f}" font-size="12">{value:.2f}</text>')
        for index, (layout, points) in enumerate(series):
            color = colors[index]
            coordinates = " ".join(f"{75+(x-1)*85:.1f},{360-320*y/maximum:.1f}"
                                   for x, y in points)
            lines.append(f'<polyline points="{coordinates}" fill="none" stroke="{color}" stroke-width="2"/>')
            for x, y in points:
                lines.append(f'<circle cx="{75+(x-1)*85:.1f}" cy="{360-320*y/maximum:.1f}" r="4" fill="{color}"/>')
            lines.append(f'<text x="{560}" y="{50+22*index}" fill="{color}" font-size="14">{layout}</text>')
        lines.append('</svg>')
        svg.write_text("\n".join(lines))
        return svg.name
    fig, ax = plt.subplots(figsize=(7, 4))
    for layout, points in series:
        ax.plot([x for x, _ in points], [y for _, y in points],
                marker="o", label=layout)
    ax.set_xticks((1, 2, 4, 8), ("1°", "0.5°", "0.25°", "0.125°"))
    ax.set_xlabel("Horizontal resolution")
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)
    ax.legend(title="DD layout")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path.name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work_dir", type=Path)
    args = parser.parse_args()
    root = args.work_dir
    with (root / "plan.csv").open(newline="") as stream:
        plan = list(csv.DictReader(stream))
    rows = []
    rank_rows = []
    for entry in plan:
        resolution = entry["resolution"]
        directory = root / f"res_{resolution.replace('.', 'p')}"
        wind_files = directory / "wind_files.csv"
        wind_gib = ""
        if wind_files.exists():
            with wind_files.open(newline="") as stream:
                wind_gib = sum(int(row["bytes"]) for row in csv.DictReader(stream)) / 2**30
        elif (directory / "wind_data").exists():
            files = list((directory / "wind_data").glob("wind_*.nc"))
            if files:
                wind_gib = sum(path.stat().st_size for path in files) / 2**30
        for layout in LAYOUTS:
            for repeat in range(1, 4):
                matches = list(directory.glob(
                    f"n*/layout{layout}/repeat{repeat}/mode1"))
                case = matches[0] if len(matches) == 1 else None
                row = {key: "" for key in FIELDS}
                row.update(resolution=resolution, nx=entry["nx"],
                           ny=entry["ny"], nz=entry["nz"],
                           layout=layout, repeat=repeat,
                           status="missing", wind_files_gib=wind_gib)
                if case:
                    result = case / "result.csv"
                    if result.exists():
                        with result.open(newline="") as stream:
                            data = next(csv.DictReader(stream))
                        row.update(status="complete", particles=data["particles"],
                            total_seconds=data["total_runtime"],
                            dd_seconds=data["dd_runtime"],
                            sent_particle_events=data["sent_particle_events"],
                            migration_fraction=data["migration_fraction"])
                    else:
                        row["status"] = "failed_or_incomplete"
                    row["input_seconds"] = input_timer(case / "trac.log")
                    io_data = traces(case / "io_trace")
                    row["io_trace_status"] = "complete" if io_data else "missing_or_incomplete"
                    row.update(io_data)
                    timing = rank_timers(case / "trac.log")
                    for rank_row in rank_traces(case / "io_trace", resolution,
                                                layout, repeat):
                        rank_row.update(timing.get(rank_row["rank"], {}))
                        rank_rows.append(rank_row)
                    row["peak_gpu_gib"] = gpu_peak(case / "gpu_memory.csv")
                    row["gpu_sample_status"] = ("complete" if row["peak_gpu_gib"] != ""
                                                else "missing")
                    if row["requested_bytes_all_ranks"] != "" and row["input_seconds"]:
                        row["requested_gib_per_input_second"] = (
                            row["requested_bytes_all_ranks"] / 2**30 /
                            row["input_seconds"])
                rows.append(row)
    with (root / "io_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    if rank_rows:
        with (root / "io_rank_summary.csv").open("w", newline="") as stream:
            fields = list(dict.fromkeys(key for row in rank_rows for key in row))
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rank_rows)
    figures = []
    for field, ylabel, filename in (
        ("input_seconds", "Input timer (s)", "input_vs_resolution.png"),
        ("total_seconds", "Total runtime (s)", "total_vs_resolution.png"),
        ("core_seconds_max_rank", "NetCDF core reads, max rank (s)", "core_vs_resolution.png"),
        ("halo_seconds_max_rank", "NetCDF periodic halo reads, max rank (s)", "halo_vs_resolution.png"),
    ):
        figure = plot(rows, field, ylabel, root / filename)
        if figure:
            figures.append(figure)
    lines = ["# GPU-DD NetCDF I/O scaling at 137 vertical levels", "",
             "The synthetic vertical coordinate spans 0–60 km. All cases use "
             "the same particle positions, six-hour simulation, four GPUs, "
             "DD exchange mode 1, and no diffusion.", "",
             "| Resolution | Wind files (GiB) | Layout | Complete / 3 | "
             "Median input (s) | Median total (s) | Median core / halo / close (s) | "
             "Requested GiB/s | Peak host / GPU (GiB) |",
             "|---|---:|---|---:|---:|---:|---:|---:|---:|"]
    for entry in plan:
        resolution = entry["resolution"]
        for layout in LAYOUTS:
            group = [row for row in rows if row["resolution"] == resolution
                     and row["layout"] == layout and row["status"] == "complete"]
            median = lambda key: (f"{statistics.median(float(row[key]) for row in group if row[key] != ''):.3f}"
                if any(row[key] != "" for row in group) else "NA")
            all_rows = [row for row in rows if row["resolution"] == resolution
                        and row["layout"] == layout]
            sizes = (f"{all_rows[0]['wind_files_gib']:.2f}"
                     if all_rows and all_rows[0]["wind_files_gib"] != "" else "NA")
            lines.append(f"| {resolution}° | {sizes} | {layout} | {len(group)} | "
                         f"{median('input_seconds')} | {median('total_seconds')} | "
                         f"{median('core_seconds_max_rank')} / "
                         f"{median('halo_seconds_max_rank')} / "
                         f"{median('close_seconds_max_rank')} | "
                         f"{median('requested_gib_per_input_second')} | "
                         f"{median('peak_host_gib_per_rank')} / "
                         f"{median('peak_gpu_gib')} |")
    lines += ["", "The NetCDF profiler records the duration of each float "
              "hyperslab call and `nc_close()` per rank. A read with a "
              "two-column longitude extent is classified as a periodic halo "
              "read. These are call times, not a decomposition of the full "
              "MPTRAC input timer, which also includes particle input and "
              "other meteorological preparation. Requested throughput is the sum of "
              "hyperslab bytes across ranks divided by the maximum input "
              "timer; file-system bytes may differ due to caching and HDF5.", "",
              "Peak host RSS is measured per rank by `getrusage`; peak GPU "
              "memory is the largest one-second `nvidia-smi` sample. Missing "
              "measurements are left blank in `io_summary.csv`.", ""]
    lines += ["## Interpretation", ""]
    for entry in plan:
        resolution = entry["resolution"]
        groups = {}
        for layout in LAYOUTS:
            group = [row for row in rows if row["resolution"] == resolution
                     and row["layout"] == layout and row["status"] == "complete"]
            if len(group) == 3 and all(row["input_seconds"] != "" for row in group):
                groups[layout] = group
        if "1x4" in groups and "4x1" in groups:
            med = lambda layout, key: statistics.median(
                float(row[key]) for row in groups[layout] if row[key] != "")
            input_ratio = med("4x1", "input_seconds") / med("1x4", "input_seconds")
            total_ratio = med("4x1", "total_seconds") / med("1x4", "total_seconds")
            send_ratio = med("4x1", "sent_particle_events") / med("1x4", "sent_particle_events")
            lines.append(
                f"At {resolution}°, `4x1` has {input_ratio:.2f}× the input "
                f"time and {total_ratio:.2f}× the total runtime of `1x4`; "
                f"its MPI send-event count is {send_ratio:.2f}× as large. "
                "The ratios describe association, not causation."
            )
            lines.append("")
    lines += [
        "The files are NetCDF4 with contiguous, uncompressed variables. "
        "In `(time, lev, lat, lon)` order, a full-longitude core read is "
        "contiguous within each latitude row; a narrow longitude slab "
        "touches many separated portions of the file. Zonal edge ranks "
        "also request two-column periodic halos. The measured core, halo, "
        "and close maxima may come from different ranks and must not be "
        "added as a decomposition of input time.", "",
        "Recommended next steps: (1) test a more contiguous core read "
        "window, measuring extra bytes against I/O time; (2) test whether "
        "periodic halo reads can reuse loaded data or be combined with "
        "compatible reads; (3) trace NetCDF/HDF5 and MPI calls around "
        "`nc_close()` to distinguish deferred work from rank waiting. "
        "No reader change is included in this baseline.", ""
    ]
    lines.extend(f"![{name}]({name})" for name in figures)
    (root / "report.md").write_text("\n".join(lines) + "\n")
    print(f"Wrote {root / 'io_summary.csv'} and {root / 'report.md'}")


if __name__ == "__main__":
    main()
