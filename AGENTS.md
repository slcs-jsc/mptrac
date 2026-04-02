# AGENTS.md

Guidance for coding agents working in this repository.

## Scope

These instructions apply to the entire repository rooted at this directory.

## Repository Map

- `src/`: primary C codebase, public header `mptrac.h`, build logic in `Makefile`, optional Fortran wrapper sources, and locally built executables.
- `tests/`: shell-driven regression and interoperability tests plus shared test data.
- `libs/`: third-party dependency build helpers and locally built libraries.
- `projects/`: examples, analysis scripts, web runner, and research workflows. Treat most subdirectories as user-facing experiments, not core library code.
- `docs/`: documentation assets.

## Working Rules

- Build from `src/`, not the repository root.
- Prefer minimal, focused changes. Avoid broad refactors unless required by the task.
- Preserve the existing C style and warning-clean build expectations. `src/Makefile` enables strict warnings and treats them as errors for `gcc` and `clang`.
- Do not edit generated build outputs in `src/` such as compiled executables or object files unless the task explicitly involves build artifacts. Edit the corresponding source files instead.
- When touching `projects/`, keep changes scoped to the relevant project; do not assume scripts there are interchangeable with the core test suite.

## Build And Test

Standard local workflow:

```bash
cd src
make
make check
```

CI currently validates these commands on Ubuntu 24.04 with optional compression and GRIB support enabled:

```bash
cd src
make ZSTD=1 ZFP=1 ECCODES=1
make check
make grib_test
```

If you change coverage-related logic, CI also uses:

```bash
cd src
make ZSTD=1 ZFP=1 ECCODES=1 COV=1
make coverage
```

## Dependencies

- Core builds expect GSL, netCDF, HDF5, and GNU Make.
- Optional CI-tested features depend on `zstd`, `zfp`, and `ecCodes`.
- Repository-local libraries are typically built under `libs/build/`, and `src/Makefile` links against that location by default.

## Testing Notes

- `make check` runs the main regression tests listed in `src/Makefile`.
- Some tests are standalone shell scripts under `tests/` such as `tests/grib_test/run.sh`, `tests/dd_test/run.sh`, and wrapper or KPP-specific checks.
- Prefer running the narrowest relevant test first, then the broader `make check` pass when your change affects shared core logic.

## Change Expectations

- Update documentation when behavior, build flags, or user workflows change.
- If a change affects command-line behavior or file formats, inspect related scripts in `tests/` and `projects/` for necessary updates.
- Keep paths and commands portable; CI uses non-interactive shell execution on Linux.

## Avoid

- Do not introduce new build systems unless explicitly requested. The repository standard is `make` via `src/Makefile`.
- Do not silently change default compiler flags or optional feature toggles without updating docs and validating affected targets.
