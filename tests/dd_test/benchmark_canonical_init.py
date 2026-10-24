#!/usr/bin/env python3
"""Set exact global-grid positions in existing DD benchmark binary inputs."""

import argparse
import array
from pathlib import Path
import struct


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("global_lon", type=int)
    parser.add_argument("global_lat", type=int)
    parser.add_argument("domains_lon", type=int)
    parser.add_argument("domains_lat", type=int)
    args = parser.parse_args()
    if args.global_lon % args.domains_lon or args.global_lat % args.domains_lat:
        parser.error("Global grid must divide evenly into DD domains")
    local_lon = args.global_lon // args.domains_lon
    local_lat = args.global_lat // args.domains_lat
    expected = local_lon * local_lat
    for rank in range(args.domains_lon * args.domains_lat):
        zonal, meridional = divmod(rank, args.domains_lat)
        path = args.data / f"data.{rank}" / "init.bin"
        with path.open("r+b") as stream:
            version, count = struct.unpack("=ii", stream.read(8))
            if version != 100 or count != expected:
                raise ValueError(f"Unexpected particle header in {path}")
            stream.seek(8 + 16 * count)
            for ix in range(local_lon):
                global_ix = zonal * local_lon + ix
                lon = 360.0 * (global_ix + 0.5) / args.global_lon
                array.array("d", [lon] * local_lat).tofile(stream)
            stream.seek(8 + 24 * count)
            lats = array.array("d", (
                90.0 - 180.0 * (meridional * local_lat + iy + 0.5)
                / args.global_lat for iy in range(local_lat)))
            for _ in range(local_lon):
                lats.tofile(stream)


if __name__ == "__main__":
    main()
