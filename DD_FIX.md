# DD Fix Notes

This note summarizes the recent fix for the `DD=1 MPI=1` path and explains how
sorting works in the current revision.

## Problem Summary

There were three separate issues:

1. `module_dd()` used `t` without receiving it as an argument.
2. `module_dd()` referenced `ctl->dd_sort_dt`, but that field does not exist.
   The correct control parameter is `ctl->sort_dt` / `SORT_DT`.
3. `dd_push()` had become implicitly OpenACC-only because it used `acc_malloc`,
   `acc_memcpy_*`, `acc_free`, and `#pragma acc` unguarded. This broke
   `DD=1 MPI=1` builds without `_OPENACC`.

The build failure was the visible symptom, but there was also a real logic bug
in the host fallback once that path was restored.

## Build-Level Fix

The compile fix was straightforward:

- `module_dd()` now takes `const double t`.
- The call site passes the current timestep `t` explicitly.
- `module_dd()` now uses `ctl->sort_dt` instead of the non-existent
  `ctl->dd_sort_dt`.

For `dd_push()` there are now two code paths:

- with `_OPENACC`: keep the device-based implementation
- without `_OPENACC`: use a host fallback

This restores support for `DD=1 MPI=1` builds that do not enable OpenACC.

## How Sorting Works Now

There are now two DD execution paths inside `module_dd()`:

- `dd_sort(...)`
- `dd_push(...)`

The choice is controlled by `SORT_DT`:

- if `SORT_DT > 0` and `fmod(t, SORT_DT) == 0`, use `dd_sort(...)`
- otherwise, use `dd_push(...)`

So sorting is optional and periodic. It is no longer required in every
timestep.

## What `dd_sort()` Does

`dd_sort()` is the "full reorder" path:

- compute sort keys for local particles
- reorder particles by box / target rank
- keep local particles at the front of the atmospheric buffer
- place outgoing particles directly behind them
- set `atm->np = nkeep`

This is mainly useful for ordering, locality, and performance.

## What `dd_push()` Does

`dd_push()` is the lightweight path:

- do not globally sort particles by box
- only copy outgoing particles behind the original local particle block
- compact the remaining local particles at the front

This is enough for correct particle exchange, even if no full sort is done in
that timestep.

## The Real Logic Bug

The subtle bug was in the interaction between `dd_push()` and
`dd_atm2particles()`.

Originally, `dd_atm2particles()` implicitly assumed that outgoing particles
always start at:

```c
ip = atm->np;
```

That assumption is valid after `dd_sort()`, because `dd_sort()` places outgoing
particles directly behind the compacted local block and then sets:

```c
atm->np = nkeep;
```

However, after restoring a host fallback for `dd_push()`, we also compacted the
local particles there. The outgoing particles were still stored behind the
*original* local block, but `atm->np` had already been changed to the compacted
local count.

That meant `dd_atm2particles()` started packing from the wrong location.

Consequences:

- wrong particles were packed for MPI exchange
- some particles disappeared
- `dd_test` no longer matched the reference output

## Actual Fix for the Logic Bug

The fix was to make the start index explicit.

`dd_atm2particles()` now takes an additional argument:

```c
const int ip0
```

This index tells the routine where the outgoing particle block starts.

`module_dd()` now sets that start index depending on the chosen DD path:

- after `dd_sort()`:
  - `ip0 = atm->np`
  - this is correct because send particles are directly behind the compacted
    local block
- after `dd_push()`:
  - `ip0 = np0`
  - where `np0` is the original local particle count before compaction

Then `dd_atm2particles()` packs:

```c
ip0 ... ip0 + npart - 1
```

instead of assuming that the send region always begins at `atm->np`.

## Why This Matters

This makes both paths correct:

- `dd_sort()` still works as before
- `dd_push()` can compact local particles without breaking MPI packing

So now:

- sorting is only an optimization / periodic maintenance step
- the plain push path is still logically correct
- the behavior is consistent for `THRUST=0` and `THRUST=1`

## Validation

The fix was verified with:

```bash
cd src
make DD=1 MPI=1 MPICC=mpicc.openmpi mptrac.o
make DD=1 MPI=1 THRUST=1 MPICC=mpicc.openmpi dd_test
```

In addition, the resulting behavior is now consistent regardless of whether the
sort backend is built with:

- `THRUST=0`
- `THRUST=1`

## Reference Data Update

`dd_test` also needed a reference update because output ATM filenames now
include seconds:

- old: `atm_YYYY_MM_DD_HH_MM.tab`
- new: `atm_YYYY_MM_DD_HH_MM_SS.tab`

After the DD fix, the remaining test mismatch was only due to this filename
change and the resulting current output ordering. The reference ATM files were
updated accordingly.
