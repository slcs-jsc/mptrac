# Contributing to MPTRAC

Thank you for contributing to MPTRAC, a Lagrangian atmospheric transport and dispersion model. Bug reports, documentation improvements, tests, and code contributions are welcome.

## Code of Conduct

Please follow the project's [Code of Conduct](CODE_OF_CONDUCT.md).

## Reporting issues

Open a [GitHub issue](https://github.com/slcs-jsc/mptrac/issues) and include a clear description, steps to reproduce (if relevant), the version or commit, input settings, error messages or logs, and the compiler, platform, and execution environment. For MPI or GPU issues, include the relevant rank, thread, accelerator, and build configuration.

## Submitting changes

1. Fork the repository, clone your fork, and create a focused branch.
2. Make the change, following the existing C style. Add or adjust relevant tests and documentation.
3. Run the applicable checks, inspect the diff, and commit with a descriptive message.
4. Push your branch and open a pull request explaining the motivation, affected workflows, tests performed, and any intended numerical or reference-data changes.

For example, from the source directory, run the default regression suite **without parallel make**:

```bash
make check
```

See [Testing and code quality](docs/manual/testing.md) for targeted tests, reference-data review, CI, GPU nightly builds, static analysis, and coverage. Do not update reference results solely to make a test pass.

## Building and coding standards

See [Building from source](docs/manual/building.md) for compilation and installation, [Code development](docs/manual/code-development.md) for formatting and documentation, and [Code structure](docs/manual/code-structure.md) for the architecture. Review the [parallelization guide](docs/manual/parallelization.md) when modifying MPI, OpenMP, or OpenACC code.

Keep changes focused and avoid unrelated formatting. Document changed behavior and public interfaces, and include relevant tests.

## License

Contributions are licensed under the [GNU General Public License v3.0](COPYING).
