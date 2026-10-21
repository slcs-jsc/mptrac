#!/bin/bash
set -e

if [[ ${MPI:-0} != 1 ]]; then
    echo "Run DD wrapper test with MPI=1" >&2
    exit 1
fi

repo_dir=$(cd ../.. && pwd)
test_root=$(mktemp -d)
trap 'rm -rf "$test_root"' EXIT

# Generate the DD fixture in an isolated tree so existing test data stays intact.
mkdir -p "$test_root/tests/dd_test" "$test_root/tests/wrapper_test/data"
ln -s "$repo_dir/src" "$test_root/src"
ln -s "$repo_dir/libs" "$test_root/libs"
cp "$repo_dir/tests/dd_test/run.sh" "$test_root/tests/dd_test/run.sh"
cp "$repo_dir/tests/wrapper_test/find-vars.sh" "$test_root/tests/wrapper_test/find-vars.sh"
export LD_LIBRARY_PATH="$repo_dir/libs/build/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1

echo "Compare C and Fortran structures..."
if ! (cd "$test_root/tests/wrapper_test" && ./find-vars.sh) \
    > "$test_root/struct.log" 2>&1; then
    tail -n 40 "$test_root/struct.log"
    exit 1
fi

echo "Build wind generator for the global grid..."
"${CC:-gcc}" -O2 -fopenmp -DNP=32 -DNQ=5 -DEX=360 -DEY=181 -DEP=90 \
    '-DVERSION="wrapper-dd-test"' \
    -I "$repo_dir/src" -I "$repo_dir/libs/build/include" \
    -o "$test_root/wind" "$repo_dir/src/wind.c" "$repo_dir/src/mptrac.c" \
    -L "$repo_dir/libs/build/lib" -lgsl -lgslcblas -lnetcdf -lm

echo "Generate DD fixture..."
if ! (cd "$test_root/tests/dd_test" && WIND_EXEC="$test_root/wind" \
    DD_EXCHANGE_MODE=1 ./run.sh skip_run skip_compare) \
    > "$test_root/fixture.log" 2>&1; then
    tail -n 40 "$test_root/fixture.log"
    if [[ -f "$test_root/tests/dd_test/data/wind_generation.log" ]]; then
        tail -n 40 "$test_root/tests/dd_test/data/wind_generation.log"
    fi
    exit 1
fi

data_dir="$test_root/tests/dd_test/data"

echo "Run C DD reference..."
if ! (cd "$data_dir" && mpirun --oversubscribe -np 9 "$repo_dir/src/trac" \
    dirlist.tab config.ctl init.nc ATM_BASENAME atm) > "$test_root/c.log" 2>&1; then
    tail -n 40 "$test_root/c.log"
    exit 1
fi

mkdir "$test_root/reference"
for subdir in "$data_dir"/data.*; do
    mkdir "$test_root/reference/$(basename "$subdir")"
    mv "$subdir"/atm_*.tab "$test_root/reference/$(basename "$subdir")/"
done

echo "Run Fortran DD wrapper..."
if ! (cd "$data_dir" && mpirun --oversubscribe -np 9 "$repo_dir/src/trac_fortran" \
    dirlist.tab config.ctl init.nc ATM_BASENAME atm) > "$test_root/fortran.log" 2>&1; then
    tail -n 40 "$test_root/fortran.log"
    exit 1
fi

echo "Compare particle tables..."
count=0
for reference in "$test_root"/reference/data.*/atm_*.tab; do
    relative=${reference#"$test_root/reference/"}
    actual="$data_dir/$relative"
    if [[ ! -f $actual ]] || ! diff -q <(sort "$reference") <(sort "$actual"); then
        echo "Different or missing: $relative" >&2
        exit 1
    fi
    count=$((count + 1))
done

for actual in "$data_dir"/data.*/atm_*.tab; do
    relative=${actual#"$data_dir/"}
    if [[ ! -f "$test_root/reference/$relative" ]]; then
        echo "Unexpected particle table: $relative" >&2
        exit 1
    fi
done

if [[ $count == 0 ]]; then
    echo "No particle tables generated" >&2
    exit 1
fi
echo "C and Fortran DD outputs match ($count particle tables)."
