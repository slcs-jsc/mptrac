#!/usr/bin/env python3
"""Summarize the 137-level DD NetCDF I/O scaling runs."""

import argparse
import csv
import html
import re
import statistics
from pathlib import Path


LAYOUTS = ("1x4", "2x2", "4x1")
DD_COMPONENTS = (
    ("dd_classification_seconds", ("dd_gpu_classify_seconds", "dd_count_outgoing_seconds")),
    ("dd_compaction_seconds", ("dd_gpu_compact_seconds",)),
    ("dd_packing_seconds", ("dd_gpu_pack_seconds", "dd_atm2particles_seconds")),
    ("dd_send_d2h_seconds", ("dd_gpu_send_d2h_seconds",)),
    ("dd_particle_communication_seconds", ("dd_communicate_particles_seconds",)),
)
RANK_TIMERS = {
    "GROUP_INPUT", "GROUP_DD", "TOTAL", "DD_GPU_CLASSIFY",
    "DD_COUNT_OUTGOING", "DD_GPU_COMPACT", "DD_GPU_PACK",
    "DD_ATM2PARTICLES", "DD_GPU_SEND_D2H", "DD_COMMUNICATE_PARTICLES",
}
FIELDS = (
    "resolution", "nx", "ny", "nz", "layout", "repeat", "read_mode", "status",
    "io_trace_status", "io_profile_source", "gpu_sample_status",
    "particles", "total_seconds", "total_seconds_mean_rank",
    "input_seconds", "input_seconds_mean_rank",
    "dd_seconds", "dd_seconds_mean_rank",
    "core_seconds_max_rank", "halo_seconds_max_rank", "close_seconds_max_rank",
    "global_seconds_max_rank", "global_seconds_mean_rank",
    "file_elapsed_seconds_max_rank", "file_elapsed_seconds_mean_rank",
    "open_seconds_max_rank", "open_seconds_mean_rank",
    "netcdf_read_seconds_max_rank", "netcdf_read_seconds_mean_rank",
    "core_seconds_mean_rank", "halo_seconds_mean_rank",
    "close_seconds_mean_rank", "empty_seconds_max_rank",
    "empty_seconds_mean_rank",
    "mpi_distribute_seconds_max_rank", "mpi_distribute_seconds_mean_rank",
    "mpi_pack_seconds_max_rank", "mpi_pack_seconds_mean_rank",
    "requested_bytes_all_ranks", "requested_gib_per_input_second",
    "core_calls_all_ranks", "halo_calls_all_ranks", "core_shape", "halo_shape",
    "global_calls_all_ranks", "global_bytes_all_ranks",
    "empty_calls_all_ranks", "netcdf_read_calls_all_ranks",
    "core_bytes_all_ranks", "halo_bytes_all_ranks",
    "peak_host_gib_per_rank", "peak_gpu_gib", "wind_files_gib",
    "sent_particle_events", "migration_fraction",
    "dd_classification_seconds_max_rank", "dd_classification_seconds_mean_rank",
    "dd_compaction_seconds_max_rank", "dd_compaction_seconds_mean_rank",
    "dd_packing_seconds_max_rank", "dd_packing_seconds_mean_rank",
    "dd_send_d2h_seconds_max_rank", "dd_send_d2h_seconds_mean_rank",
    "dd_particle_communication_seconds_max_rank",
    "dd_particle_communication_seconds_mean_rank",
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
        match = re.match(r"\s*(\d+): TIMER_([A-Z0-9_]+) =\s*([\d.eE+-]+) s", line)
        if match:
            rank, name, value = match.groups()
            if name not in RANK_TIMERS:
                continue
            detail = re.search(r"mean=\s*([\d.eE+-]+) s, max=.*?n=\s*(\d+)\)", line)
            if detail:
                value = float(detail.group(1)) * int(detail.group(2))
            timers.setdefault(int(rank), {})[name.lower() + "_seconds"] = float(value)
    return timers


def trace_records(path, log):
    """Return per-rank, per-file records from native logs or legacy CSV."""
    records = {}
    if log.exists():
        for line in log.read_text(errors="replace").splitlines():
            match = re.search(r"\bMET_IO rank=(\d+) file=(.*?) (open_s=.*)", line)
            if not match:
                continue
            rank = int(match.group(1))
            values = dict(item.split("=", 1) for item in match.group(3).split())
            values["file"] = match.group(2)
            values.setdefault("global_s", "0")
            values.setdefault("global_calls", "0")
            values.setdefault("global_bytes", "0")
            for current, previous in (("met_nc_read_s", "read_s"),
                                      ("met_mpi_distribute_s", "mpi_distribute_s"),
                                      ("met_mpi_pack_s", "mpi_pack_s")):
                if current not in values and previous in values:
                    values[current] = values.pop(previous)
            for source, target in (("open_s", "open_seconds"),
                                   ("met_nc_read_s", "netcdf_read_seconds"),
                                   ("core_s", "core_seconds"),
                                   ("halo_s", "halo_seconds"),
                                   ("global_s", "global_seconds"),
                                   ("file_elapsed_s", "file_elapsed_seconds"),
                                   ("empty_s", "empty_seconds"),
                                   ("close_s", "close_seconds"),
                                   ("met_mpi_distribute_s", "mpi_distribute_seconds"),
                                   ("met_mpi_pack_s", "mpi_pack_seconds")):
                if source in values:
                    values[target] = values.pop(source)
            for category in ("core", "halo"):
                shape = values.pop(f"{category}_shape").split("x")
                for dim, size in zip(("lev", "lat", "lon"), shape):
                    values[f"{category}_{dim}"] = size
            records.setdefault(rank, []).append(values)
    if records:
        return records, "native"
    for file in sorted(path.glob("rank_*.csv")):
        with file.open(newline="") as stream:
            items = list(csv.DictReader(stream))
        for item in items:
            if int(item["halo_bytes"]) == 0:
                item["empty_seconds"] = str(float(item.get("empty_seconds", 0)) +
                                            float(item["halo_seconds"]))
                item["empty_calls"] = str(int(item.get("empty_calls", 0)) +
                                          int(item["halo_calls"]))
                item["halo_seconds"] = "0"
                item["halo_calls"] = "0"
            else:
                item.setdefault("empty_seconds", "0")
                item.setdefault("empty_calls", "0")
            item["netcdf_read_seconds"] = str(sum(float(item[key]) for key in
                ("core_seconds", "halo_seconds", "empty_seconds")))
            item["global_seconds"] = "0"
            item["global_calls"] = "0"
            item["global_bytes"] = "0"
        records[int(file.stem.split("_")[-1])] = items
    return records, "legacy"


def traces(records, ranks=4, files_per_rank=7):
    if set(records) != set(range(ranks)) or any(
            len(items) != files_per_rank for items in records.values()):
        return {}
    per_rank = [records[rank] for rank in range(ranks)]
    def sums(key):
        if any(key not in record for items in per_rank for record in items):
            return None
        return [sum(float(record[key]) for record in items) for items in per_rank]
    def aggregate(key):
        values = sums(key)
        return ((max(values), statistics.mean(values)) if values is not None
                else ("", ""))
    shapes = lambda prefix: ";".join(sorted({"x".join(
        record[f"{prefix}_{dim}"] for dim in ("lev", "lat", "lon"))
        for items in per_rank for record in items
        if record[f"{prefix}_lev"] != "0"}))
    result = {}
    for source, target in (("open_seconds", "open_seconds"),
                           ("netcdf_read_seconds", "netcdf_read_seconds"),
                           ("core_seconds", "core_seconds"),
                           ("halo_seconds", "halo_seconds"),
                           ("global_seconds", "global_seconds"),
                           ("file_elapsed_seconds", "file_elapsed_seconds"),
                           ("empty_seconds", "empty_seconds"),
                           ("close_seconds", "close_seconds"),
                           ("mpi_distribute_seconds", "mpi_distribute_seconds"),
                           ("mpi_pack_seconds", "mpi_pack_seconds")):
        maximum, mean = aggregate(source)
        result[f"{target}_max_rank"] = maximum
        result[f"{target}_mean_rank"] = mean
    for source, target in (("core_calls", "core_calls_all_ranks"),
                           ("halo_calls", "halo_calls_all_ranks"),
                           ("global_calls", "global_calls_all_ranks"),
                           ("empty_calls", "empty_calls_all_ranks"),
                           ("core_bytes", "core_bytes_all_ranks"),
                           ("halo_bytes", "halo_bytes_all_ranks"),
                           ("global_bytes", "global_bytes_all_ranks")):
        values = sums(source)
        result[target] = int(sum(values)) if values is not None else ""
    calls = (result["core_calls_all_ranks"], result["halo_calls_all_ranks"],
             result["global_calls_all_ranks"], result["empty_calls_all_ranks"])
    result["netcdf_read_calls_all_ranks"] = (sum(calls) if all(
        value != "" for value in calls) else "")
    result["requested_bytes_all_ranks"] = (
        result["core_bytes_all_ranks"] + result["halo_bytes_all_ranks"]
        + result["global_bytes_all_ranks"])
    result["core_shape"] = shapes("core")
    result["halo_shape"] = shapes("halo")
    result["peak_host_gib_per_rank"] = max(float(record["peak_rss_kib"])
        for items in per_rank for record in items) / 2**20
    return result


def rank_traces(records, resolution, layout, repeat, files_per_rank=7):
    rows = []
    for rank, items in sorted(records.items()):
        if len(items) != files_per_rank:
            continue
        row = {"resolution": resolution, "layout": layout, "repeat": repeat,
               "rank": rank}
        for field in ("open_seconds", "netcdf_read_seconds", "core_seconds",
                      "halo_seconds", "global_seconds", "empty_seconds",
                      "close_seconds", "file_elapsed_seconds",
                      "mpi_distribute_seconds", "mpi_pack_seconds",
                      "core_calls", "halo_calls", "global_calls", "empty_calls",
                      "core_bytes", "halo_bytes", "global_bytes"):
            if all(field in record for record in items):
                convert = int if field.endswith(("_calls", "_bytes")) else float
                row[field] = sum(convert(record[field]) for record in items)
        if all(field in row for field in ("core_calls", "halo_calls",
                                          "global_calls", "empty_calls")):
            row["netcdf_read_calls"] = sum(row[field] for field in
                ("core_calls", "halo_calls", "global_calls", "empty_calls"))
        row["requested_bytes"] = (row["core_bytes"] + row["halo_bytes"]
                                  + row["global_bytes"])
        row["peak_host_gib"] = max(float(record["peak_rss_kib"])
                                    for record in items) / 2**20
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
    parser.add_argument("--files-per-run", type=int, default=7)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--gpus", type=int, default=4)
    parser.add_argument("--read-mode", type=int, choices=(0, 1), default=0)
    args = parser.parse_args()
    if min(args.files_per_run, args.repeats, args.ranks) < 1:
        parser.error("file count, repeats, and ranks must be positive")
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
            for repeat in range(1, args.repeats + 1):
                matches = list(directory.glob(
                    f"n*/layout{layout}/repeat{repeat}/mode1"))
                case = matches[0] if len(matches) == 1 else None
                row = {key: "" for key in FIELDS}
                row.update(resolution=resolution, nx=entry["nx"],
                           ny=entry["ny"], nz=entry["nz"],
                           layout=layout, repeat=repeat, read_mode=args.read_mode,
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
                    records, source = trace_records(case / "io_trace",
                                                    case / "trac.log")
                    io_data = traces(records, args.ranks, args.files_per_run)
                    row["io_trace_status"] = "complete" if io_data else "missing_or_incomplete"
                    row["io_profile_source"] = source if io_data else ""
                    row.update(io_data)
                    timing = rank_timers(case / "trac.log")
                    if set(timing) == set(range(args.ranks)):
                        for source, target in (("total_seconds", "total_seconds_mean_rank"),
                                               ("group_input_seconds", "input_seconds_mean_rank"),
                                               ("group_dd_seconds", "dd_seconds_mean_rank")):
                            values = [timing[rank].get(source)
                                      for rank in range(args.ranks)]
                            if all(value is not None for value in values):
                                row[target] = statistics.mean(values)
                        for target, sources in DD_COMPONENTS:
                            for source in sources:
                                values = [timing[rank].get(source)
                                          for rank in range(args.ranks)]
                                if all(value is not None for value in values):
                                    row[target + "_max_rank"] = max(values)
                                    row[target + "_mean_rank"] = statistics.mean(values)
                                    break
                    profiles = {item["rank"]: item for item in rank_traces(
                        records, resolution, layout, repeat, args.files_per_run)}
                    for rank in sorted(set(profiles) | set(timing)):
                        rank_row = profiles.get(rank, {
                            "resolution": resolution, "layout": layout,
                            "repeat": repeat, "rank": rank})
                        rank_row.update(timing.get(rank, {}))
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
        ("global_seconds_max_rank", "NetCDF global 2-D reads, max rank (s)", "global_vs_resolution.png"),
    ):
        figure = plot(rows, field, ylabel, root / filename)
        if figure:
            figures.append(figure)
    platform = "GPU-DD" if args.gpus else "CPU-DD"
    lines = [f"# {platform} NetCDF I/O scaling at 137 vertical levels", "",
             "The synthetic vertical coordinate spans 0–60 km. All cases use "
             f"the same particle positions, {args.files_per_run - 1}-hour "
             f"simulation, {args.ranks} MPI ranks, {args.gpus} GPUs, "
             f"DD exchange mode 1, NetCDF read mode {args.read_mode}, "
             "and no diffusion.", "",
             f"| Resolution | Wind files (GiB) | Layout | Complete / {args.repeats} | "
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
    lines += ["", "## Call-level measurements", "",
              f"Each entry is the median of {args.repeats} run-level summaries. Time "
              "pairs are maximum / mean of the **rank-local accumulated call "
              "times**; they are not elapsed times for a collective phase. "
              "A missing measurement is `NA`.", "",
              "| Grid | Layout | Open max / mean (s) | NetCDF read max / mean (s) | "
              "MPI redistribute max / mean (s) | MPI pack max / mean (s) | "
              "Core / halo / global / empty calls | "
              "Core / halo / global requested GiB |",
              "|---|---|---:|---:|---:|---:|---:|---:|"]
    for entry in plan:
        resolution = entry["resolution"]
        for layout in LAYOUTS:
            group = [row for row in rows if row["resolution"] == resolution
                     and row["layout"] == layout and row["status"] == "complete"]
            def median(key, scale=1):
                values = [float(row[key]) / scale for row in group
                          if row[key] != ""]
                return f"{statistics.median(values):.3f}" if values else "NA"
            def pair(prefix):
                return (f"{median(prefix + '_max_rank')} / "
                        f"{median(prefix + '_mean_rank')}")
            lines.append(
                f"| {resolution}° | {layout} | {pair('open_seconds')} | "
                f"{pair('netcdf_read_seconds')} | "
                f"{pair('mpi_distribute_seconds')} | "
                f"{pair('mpi_pack_seconds')} | "
                f"{median('core_calls_all_ranks')} / "
                f"{median('halo_calls_all_ranks')} / "
                f"{median('global_calls_all_ranks')} / "
                f"{median('empty_calls_all_ranks')} | "
                f"{median('core_bytes_all_ranks', 2**30)} / "
                f"{median('halo_bytes_all_ranks', 2**30)} / "
                f"{median('global_bytes_all_ranks', 2**30)} |")
    lines += ["", "## DD and I/O timing breakdown", "",
              "Values are seconds. Call-level and DD phase timers are summed "
              "over the run on each rank; the maximum rank-local value is "
              "then taken for each component, followed by the median across "
              "completed repetitions. Group and total timers use their "
              "direct MPTRAC values. The "
              "meteorological file phase is `MET_IO file_elapsed_s`; older "
              "traces without it show N/A. The MPTRAC input group includes "
              "particle input and other input work. Component maxima may "
              "come from different ranks and must not be added.", ""]
    components = (
        ("NetCDF open", "open_seconds_max_rank"),
        ("NetCDF reads", "netcdf_read_seconds_max_rank"),
        ("NetCDF close", "close_seconds_max_rank"),
        ("MPI meteorological redistribution", "mpi_distribute_seconds_max_rank"),
        ("Meteorological packing", "mpi_pack_seconds_max_rank"),
        ("Total meteorological input (file phase)", "file_elapsed_seconds_max_rank"),
        ("MPTRAC input group", "input_seconds"),
        ("DD particle classification", "dd_classification_seconds_max_rank"),
        ("DD particle compaction", "dd_compaction_seconds_max_rank"),
        ("DD particle packing", "dd_packing_seconds_max_rank"),
        ("DD outgoing device-to-host", "dd_send_d2h_seconds_max_rank"),
        ("DD particle communication phase", "dd_particle_communication_seconds_max_rank"),
        ("Total particle DD", "dd_seconds"),
        ("Total simulation", "total_seconds"),
    )
    for entry in plan:
        resolution = entry["resolution"]
        lines += [f"### {resolution}° — NetCDF read mode {args.read_mode}", "",
                  "| Component | 1x4 | 2x2 | 4x1 |",
                  "|---|---:|---:|---:|"]
        for label, key in components:
            cells = []
            for layout in LAYOUTS:
                values = [float(row[key]) for row in rows
                          if row["resolution"] == resolution
                          and row["layout"] == layout
                          and row["status"] == "complete" and row[key] != ""]
                cells.append(f"{statistics.median(values):.3f}"
                             if values else "N/A")
            lines.append(f"| {label} | {' | '.join(cells)} |")
        lines.append("")
    lines += ["`DD_GPU_CLASSIFY` (or CPU `DD_COUNT_OUTGOING`) supplies "
              "classification; `DD_GPU_COMPACT`, `DD_GPU_PACK`, and "
              "`DD_GPU_SEND_D2H` supply the corresponding GPU rows. "
              "`DD_ATM2PARTICLES` is the packing fallback for the sort/push "
              "path. `DD_COMMUNICATE_PARTICLES` includes host preparation "
              "as well as MPI calls. `GROUP_DD` and `TOTAL` come directly "
              "from MPTRAC timers. N/A means the measurement is absent, "
              "not zero; rank-local timer values remain in "
              "`io_rank_summary.csv`.", ""]
    lines += ["", "The optional native DD profiler measures `nc_open_par()`, "
              "each float core/periodic-halo/global-slice read, "
              "and `nc_close()` with a local monotonic wall clock. "
              "Zero-element collective read calls have their own `empty` "
              "category. `MET_NC_READ` is the sum of core, halo, global, "
              "and empty NetCDF call durations on a rank. MPI redistribution "
              "and CPU pack counters are nonzero only in read mode 1. "
              "`file_elapsed_s` spans each rank's complete NetCDF file read "
              "and redistribution without a synchronization barrier. "
              "Archived LD_PRELOAD traces lack open and MPI timings. The "
              "current interposer records open time but still classifies "
              "halos by the last hyperslab dimension. Missing fields "
              "remain `NA`.", "",
              "The MPTRAC input timer also includes particle input and "
              "meteorological preparation. It cannot be decomposed by "
              "adding independent rank maxima. Requested throughput is the "
              "sum of hyperslab bytes across ranks divided by the maximum "
              "input timer; actual filesystem bytes may differ due to "
              "caching and HDF5. Per-rank values are in `io_rank_summary.csv`.", "",
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
            if len(group) == args.repeats and all(
                    row["input_seconds"] != "" for row in group):
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
    if args.read_mode == 0:
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
    else:
        lines += [
            "The global 2-D reader requests complete horizontal slices and "
            "uses MPI to distribute the local DD windows. Compare its "
            "global NetCDF reads and MPI redistribution using the measured "
            "rank-local values above; neither component alone is the total "
            "meteorological input time.", ""
        ]
    lines.extend(f"![{name}]({name})" for name in figures)
    (root / "report.md").write_text("\n".join(lines) + "\n")
    print(f"Wrote {root / 'io_summary.csv'} and {root / 'report.md'}")


if __name__ == "__main__":
    main()
