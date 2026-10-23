#!/bin/bash
set -euo pipefail

# Run on a GPU node with DD_USE_GPU=1 DD_COMPILER=nvc DD_MPICC=mpicc.

script_dir=$(cd "$(dirname "$0")" && pwd)
mptrac_dir=$(cd "$script_dir/../.." && pwd)
work_dir=${DD_UVWP_WORK_DIR:-$(mktemp -d /tmp/mptrac_dd_uvwp.XXXXXX)}
mkdir -p "$work_dir"
echo "[INFO] uvwp test directory: $work_dir"

mpi_cc=${DD_MPICC:-/usr/bin/mpicc}
mpi_run=${DD_MPIRUN:-/usr/bin/mpirun}
compiler=${DD_COMPILER:-gcc}
gpu=${DD_USE_GPU:-0}
defs='-DNP=32 -DEX=8 -DEY=8 -DEP=4 -DNQ=3'

if [[ $compiler == nvc ]]; then
    cflags=(-mp)
    if [[ $gpu -eq 1 ]]; then cflags+=(-acc -gpu=cc80,cc90); fi
else
    cflags=(-fopenmp)
    if [[ $gpu -eq 1 ]]; then cflags+=(-fopenacc -foffload=-lm); fi
fi

make -C "$mptrac_dir/src" clean > "$work_dir/build.log" 2>&1
make -C "$mptrac_dir/src" mptrac.o DD=1 MPI=1 GPU="$gpu" \
    COMPILER="$compiler" MPICC="$mpi_cc" STATIC=0 DEFINES="$defs" \
    >> "$work_dir/build.log" 2>&1
"$mpi_cc" -O2 $defs -DMPI -DDD "${cflags[@]}" \
    -I "$mptrac_dir/src" -I "$mptrac_dir/libs/build/include" \
    "$script_dir/uvwp.c" "$mptrac_dir/src/mptrac.o" \
    -L "$mptrac_dir/libs/build/lib" -lgsl -lgslcblas -lnetcdf -lm \
    -o "$work_dir/uvwp" >> "$work_dir/build.log" 2>&1

if [[ ${DD_UVWP_COMPILE_ONLY:-0} == 1 ]]; then
    echo "[INFO] uvwp MPI/OpenACC test compiled"
    exit 0
fi

export LD_LIBRARY_PATH="$mptrac_dir/libs/build/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1
for path in sort push compact; do
    "$mpi_run" --oversubscribe -np 2 "$work_dir/uvwp" "$path" \
        > "$work_dir/$path.log" 2>&1
    echo "[INFO] uvwp migration passed: $path"
done
