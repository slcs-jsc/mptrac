#!/usr/bin/env python3
"""Split a global MPTRAC binary init file into per-rank files by DD ownership.

Usage: init_partition.py init_global.bin out_dir nx ny [layouts]

layouts is a comma-separated list such as "1x4,2x2,4x1" (default).
Writes out_dir/init_layout<L>/data/data.<r>/init.bin for each layout and
prints per-layout per-rank counts as "layout c0,c1,c2,c3".

Ownership replicates dd_assign_subdomains(): rank -> zonal = rank // dlat,
merid = rank % dlat; index blocks split lon/lat grid arrays with the last
rank absorbing the remainder; a particle belongs to the rank owning
locate_reg(lon) x locate_irr(lat) on the met file's regular coordinate grids.
"""
import array
import pathlib
import struct
import sys

NQ = 5
COLS = 4 + NQ


def locate_reg(step, x0, n, x):
    i = int((x - x0) / step)
    if i < 0:
        return 0
    if i > n - 2:
        return n - 2
    return i


def locate_irr(xx, n, x):
    ilo, ihi = 0, n - 1
    while ihi > ilo + 1:
        i = (ihi + ilo) >> 1
        if xx[i] <= x:
            ihi = i
        else:
            ilo = i
    return ilo


def grids(nx, ny):
    lon = [360.0 * i / float(nx) for i in range(nx)]
    lat = [90.0 - 180.0 * j / float(ny - 1) for j in range(ny)]
    return lon, lat


def read_columns(path, count):
    with path.open("rb") as stream:
        version, npart = _header(stream, path)
        if npart != count:
            raise SystemExit(f"[ERROR] {path}: header count {npart} != {count}")
        columns = []
        for _ in range(COLS):
            column = array.array("d")
            column.fromfile(stream, count)
            columns.append(column)
        final = array.array("i")
        final.fromfile(stream, 1)
        if len(final) != 1 or final[0] != 999:
            raise SystemExit(f"[ERROR] {path}: missing final flag")
    return columns


def _header(stream, path):
    version, count = struct.unpack("=ii", stream.read(8))
    if version != 100:
        raise SystemExit(f"[ERROR] {path}: unexpected version {version}")
    return version, count


def write_columns(path, count, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as stream:
        stream.write(struct.pack("=ii", 100, count))
        for column in columns:
            column.tofile(stream)
        stream.write(struct.pack("=i", 999))


def main():
    src = pathlib.Path(sys.argv[1])
    out_root = pathlib.Path(sys.argv[2])
    nx, ny = int(sys.argv[3]), int(sys.argv[4])
    layouts = sys.argv[5] if len(sys.argv) > 5 else "1x4,2x2,4x1"
    lon_step = 360.0 / float(nx)
    lat_grid = grids(nx, ny)[1]
    lat_cache = {}
    with src.open("rb") as stream:
        _, count = _header(stream, src)
    columns = read_columns(src, count)
    ids = columns[4]
    if sorted(ids) != list(range(count)):
        raise SystemExit("[ERROR] global IDs are not exactly 0..count-1")
    lon, lat = columns[2], columns[3]
    norm = []
    for i in range(count):
        lo = lon[i] % 360.0
        la = lat[i]
        if la > 90.0:
            la = 180.0 - la
            lo = (lo + 180.0) % 360.0
        elif la < -90.0:
            la = -180.0 - la
            lo = (lo + 180.0) % 360.0
        norm.append((lo, la))
    ix = [locate_reg(lon_step, 0.0, nx, norm[i][0]) for i in range(count)]
    iy = []
    for i in range(count):
        v = norm[i][1]
        j = lat_cache.get(v)
        if j is None:
            j = locate_irr(lat_grid, ny, v)
            lat_cache[v] = j
        iy.append(j)
    for layout in layouts.split(","):
        dlon, dlat = (int(part) for part in layout.lower().split("x"))
        nblk_x = nx // dlon
        nblk_y = ny // dlat
        rix = [min(i // nblk_x, dlon - 1) for i in ix]
        riy = [min(i // nblk_y, dlat - 1) for i in iy]
        groups = [[] for _ in range(dlon * dlat)]
        for i in range(count):
            groups[rix[i] * dlat + riy[i]].append(i)
        counts = []
        for rank in range(dlon * dlat):
            keep = groups[rank]
            sub = [array.array("d", (col[i] for i in keep)) for col in columns]
            write_columns(out_root / f"init_layout{layout}" / "data" /
                          f"data.{rank}" / "init.bin", len(keep), sub)
            counts.append(len(keep))
        if sum(counts) != count:
            raise SystemExit(f"[ERROR] layout {layout}: partition lost particles")
        print(layout, ",".join(map(str, counts)))


if __name__ == "__main__":
    main()
