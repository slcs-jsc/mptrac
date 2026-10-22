#!/bin/bash
set -euo pipefail

# Compare complete DD snapshots by particle ID (column 5), independently of
# OpenACC's atomic compaction order. Call from tests/dd_test.
actual=${1:-data}
reference=${2:-data.ref}
if [[ ! -d "$actual" || ! -d "$reference" ]]; then
    echo "[ERROR] Missing DD output directory: $actual or $reference" >&2
    exit 1
fi

status=0
while IFS= read -r -d '' ref_file; do
    rel_path=${ref_file#"$reference"/}
    data_file=$actual/$rel_path
    if [[ ! -f "$data_file" ]]; then
        echo "[ERROR] Missing file: $data_file" >&2
        status=1
    elif [[ $(basename "$ref_file") == atm_*.tab ]]; then
        if ! diff -q \
            <(awk 'NF && $1 !~ /^#/ {print}' "$data_file" | LC_ALL=C sort -k5,5n) \
            <(awk 'NF && $1 !~ /^#/ {print}' "$ref_file" | LC_ALL=C sort -k5,5n) \
            > /dev/null; then
            echo "[ERROR] Particle IDs or values differ: $data_file" >&2
            status=1
        fi
    elif ! diff -q "$data_file" "$ref_file" > /dev/null; then
        echo "[ERROR] Files differ: $data_file" >&2
        status=1
    fi
done < <(find "$reference" -type f -name '*.tab' -print0)

while IFS= read -r -d '' data_file; do
    rel_path=${data_file#"$actual"/}
    if [[ ! -f "$reference/$rel_path" ]]; then
        echo "[ERROR] Extra file: $data_file" >&2
        status=1
    fi
done < <(find "$actual" -type f -name '*.tab' -print0)

# A particle ID must occur once across all ranks at each output time.
while IFS= read -r stamp; do
    for directory in "$actual" "$reference"; do
        if ! awk '
            NF && $1 !~ /^#/ {
                if (seen[$5]++) {
                    printf "[ERROR] Duplicate particle ID %s at time %s\n", $5, $1 > "/dev/stderr"
                    bad = 1
                }
            }
            END {exit bad}
        ' "$directory"/data.*/"atm_$stamp.tab"; then
            status=1
        fi
    done
done < <(find "$reference" -type f -name 'atm_*.tab' -printf '%f\n' | sed 's/^atm_//; s/\.tab$//' | sort -u)

if [[ $status -eq 0 ]]; then
    echo "[INFO] DD snapshots match by particle ID: $actual vs $reference"
fi
exit "$status"
