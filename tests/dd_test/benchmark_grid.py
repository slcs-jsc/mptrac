#!/usr/bin/env python3
"""Print an exact DD benchmark particle grid as '<lon> <lat>'.

The grid satisfies lon * lat == target exactly, with lon a multiple of
step_lon and lat a multiple of step_lat so that every DD layout divides
the lattice evenly.  The aspect ratio stays near the historical
lat ~ sqrt(target / 2) heuristic as closely as the divisors allow.
"""

import sys


def grid(target, step_lon, step_lat):
    best = None
    lat = step_lat
    while lat * step_lat <= target * step_lat and lat <= target:
        if target % lat == 0:
            lon = target // lat
            if lon % step_lon == 0:
                error = abs(lat * 2 - target ** 0.5 * (2 ** 0.5))
                if best is None or error < best[0]:
                    best = (error, lon, lat)
        lat += step_lat
        if lat > 2 * (target ** 0.5) + step_lat:
            break
    if best is None:
        raise SystemExit(
            f"[ERROR] No exact grid for {target} particles with strides "
            f"{step_lon}x{step_lat}; choose a count divisible by "
            f"{step_lon * step_lat} with suitable divisors")
    return best[1], best[2]


def main():
    if len(sys.argv) != 4:
        raise SystemExit("Usage: benchmark_grid.py <target> <step_lon> <step_lat>")
    target, step_lon, step_lat = map(int, sys.argv[1:4])
    if min(target, step_lon, step_lat) < 1:
        raise SystemExit("[ERROR] Arguments must be positive integers")
    lon, lat = grid(target, step_lon, step_lat)
    print(lon, lat)


if __name__ == "__main__":
    main()
