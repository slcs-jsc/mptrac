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
- When changing compile-time options, run `make clean` before rebuilding; Make does not track changes to command-line flags.
- Prefer minimal, focused changes. Avoid broad refactors unless required by the task.
- Before changing behavior, check the existing implementation, callers, tests, and documented expectations. Challenge assumptions and explain non-obvious choices in the change summary.
- Keep solutions as simple as the requirement allows. Add abstractions, dependencies, or compatibility paths only when a concrete use case needs them.
- Avoid new `static` helper functions. Use a macro for a small operation when its arguments can be evaluated safely, or use a regular function declared and documented with Doxygen in the relevant header (usually `src/mptrac.h`).
- Reuse existing data structures where possible. Introduce a new `struct` only when the required state does not fit an existing one cleanly; first check the established OpenACC host/device transfers and other code paths that depend on those structures.
- When changing public declarations or data layouts in `src/mptrac.h`, check affected C callers, the Fortran wrapper in `src/mptrac_fortran.f90`, and MPI/OpenACC data paths as applicable. Run `make wrapper_test` for changes that affect the C/Fortran interface.
- Preserve the existing C style and warning-clean build expectations. `src/Makefile` enables strict warnings and treats them as errors for `gcc` and `clang`.
- Format C source and headers from `src/` with `make indent` after a successful build. This target rewrites all `*.c` and `*.h` files there, so review the diff and keep only changes relevant to the task.
- Do not edit generated build outputs in `src/` such as compiled executables or object files unless the task explicitly involves build artifacts. Edit the corresponding source files instead.
- When touching `projects/`, keep changes scoped to the relevant project; do not assume scripts there are interchangeable with the core test suite.

## Build And Test

Standard local workflow:

```bash
cd src
make
make check
```

Run `make check` without `-j`; parallel test execution can hide failures.

The main CI test workflow on Ubuntu 24.04 builds with optional compression and GRIB support, then checks the Fortran wrapper:

```bash
cd src
make ZSTD=1 LZ4=1 ZFP=1 ECCODES=1
make check
make grib_test
make wrapper ZSTD=1 LZ4=1 ZFP=1 ECCODES=1
make wrapper_test ZSTD=1 LZ4=1 ZFP=1 ECCODES=1
```

CI also runs a coverage workflow. To include optional features in a coverage report, pass their flags to `make coverage`, which cleans and rebuilds:

```bash
cd src
make coverage ZSTD=1 LZ4=1 ZFP=1 ECCODES=1
```

## Dependencies

- Core builds expect GSL, netCDF, HDF5, and GNU Make.
- Optional CI-tested features depend on `zstd`, `lz4`, `zfp`, and `ecCodes`. The Fortran wrapper also needs a Fortran compiler.
- Repository-local libraries are typically built under `libs/build/`, and `src/Makefile` links against that location by default.

## Testing Notes

- `make check` runs the main regression tests listed in `src/Makefile`. Individual tests can be run from `src/` with targets such as `make trac_test`.
- Additional targets include `grib_test`, `wrapper_test`, `dd_test`, and `kpp_test`; their scripts are under `tests/`.
- Prefer running the narrowest relevant test first, then the broader `make check` pass when your change affects shared core logic.
- Investigate test failures against the intended behavior. Update reference data only when the new results are understood and verified as correct.

## Change Expectations

- Keep the MkDocs user manual in `docs/manual/` and the Doxygen API documentation in sync with relevant code, build, and workflow changes. Build the affected documentation from `src/` with `make mkdocs` and/or `make doxygen` when the required tools are available.
- For new or revised Doxygen API comments, follow the established `/** ... */` blocks in `src/mptrac.h`: use `@brief` plus a substantive description, document parameters, return values, and behavior where relevant, and always include `@author` with the responsible author's name. Do not invent attribution. Preserve the surrounding style for other comments.
- If a change affects command-line behavior or file formats, inspect related scripts in `tests/` and `projects/` for necessary updates.
- Keep paths and commands portable; CI uses non-interactive shell execution on Linux.

## Avoid

- Do not introduce new build systems unless explicitly requested. The repository standard is `make` via `src/Makefile`.
- Do not silently change default compiler flags or optional feature toggles without updating docs and validating affected targets.
