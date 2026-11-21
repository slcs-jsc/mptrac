# Code development

Please follow these guidelines when developing new code for MPTRAC.
For issue reports and pull requests, see the [contribution guidelines](https://github.com/slcs-jsc/mptrac/blob/master/CONTRIBUTING.md).

## Code formatting

On Linux, use the `indent` tool to apply a selected formatting style
to the source code:

```
    make indent
```

Only use `indent` if the code has been compiled correctly. You may
need to run the `indent` command 2-3 times to achieve the desired
formatting.


## Documentation

Please update the [README file](https://github.com/slcs-jsc/mptrac/blob/master/README.md) as needed.

To update the [user manual](https://slcs-jsc.github.io/mptrac), use:

```
    make mkdocs
```

To update the [Doxygen manual](https://slcs-jsc.github.io/mptrac/doxygen), use:

```
    make doxygen
```

## Further reading

- [Git Handbook](https://guides.github.com/introduction/git-handbook)

- [GitHub Git Cheat Sheet](https://training.github.com/downloads/github-git-cheat-sheet.pdf)

## Related developer pages

- [Code structure](code-structure.md): architecture and major data structures.
- [Building from source](building.md): compilation and installation.
- [Testing and code quality](testing.md): regression tests, static analysis, coverage, and automated checks.
- [Parallelization](parallelization.md): MPI, OpenMP, and OpenACC.
- [Profiling and benchmarking](profiling.md): timers and CPU/GPU profiling.
- [Release management](releases.md): release policy and checklist.
- [High-level interface](high-level-interface.md) and [Fortran wrapper](fortran-wrapper.md): integration APIs.
