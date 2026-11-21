# Code development

Please follow these guidelines when developing new code for MPTRAC.

## Compilation

To remove old binaries, backup files, etc., and to clean up the source code directory:

```
    make clean
```

For faster compilation, you can use parallel compilation with `make`:

```
    make -j
```

To compile a specific binary (e.g., `trac`):

```
    make trac
```

You can change Makefile variables at the command line:

```
    make DEFINES="-DNP=100000" GPU=1 STATIC=0
```

## Testing

Always run the full test suite to verify your code changes:

```
    make check
```

Do not run tests in parallel (i.e., avoid using the `-j` option with
`make check`), as this may cause test failures to be missed.

You can run only selected tests, e.g.:

```
    make trac_test
```

If a test fails, please carefully compare the test results with the
reference data. The Linux tool `xxdiff` can be used to compare test
data and reference data number by number. The graphing utility
`gnuplot` can be used to visualize any differences. The reference data
of a test should only be updated if the new test results are
considered correct.

## Static code analysis

You can use static code analysis to automatically detect potential
problems in the code:

```
    make cppcheck
    make lizard
```

## Coverage testing

You can perform a coverage analysis to determine which parts of the
code are covered by tests:

```
    make COV=1
    make coverage
```

If you find that parts of the code are not covered by the coverage
analysis, please consider adding a new test.

## Automated testing and code analysis

After committing revised code to the GitHub repository, please check
the [GitHub Actions page](https://github.com/slcs-jsc/mptrac/actions)
to see if the automated tests were successfully passed.

Please also check the
[Codacy](https://app.codacy.com/gh/slcs-jsc/mptrac?utm_source=github.com&utm_medium=referral&utm_content=slcs-jsc/mptrac&utm_campaign=Badge_Grade_Settings)
and [Codecov](https://codecov.io/gh/slcs-jsc/mptrac) websites for test
results.

Please check the
[nightly build website](https://datapub.fz-juelich.de/slcs/mptrac/nightly_builds)
to see if the tests passed on the Juelich supercomputers, especially for
the GPU code.

Users at the Jülich Supercomputing Centre can also check the test
results at the [ESM Buildbot](https://esm-buildbot.fz-juelich.de).

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

## Installation

To copy the executables from the source directory to the `DESTDIR`
directory (default `../bin/`):

```
    make install
```

To remove the executables from the installation directory:

```
    make uninstall
```

To package the current repository state into a zip file, including the
compiled binaries, run:

```
    make dist
```

## Releases

MPTRAC and JURASSIC follow a coordinated, lightweight release process.
Both projects aim for approximately two regular releases per year,
preferably in January and July, when development activity is relatively low.
Release dates are flexible; the projects may also publish additional
releases after major development or validation milestones, such as a
comprehensive code review. The two projects should normally be released
around the same time, but a release of one project does not depend on the
readiness of the other.

Version numbers reflect the changes made, not the release date. Existing
Git tags use the historical `vX.Y` convention. These two-component tags
are retained for compatibility with previous releases; they are not full
three-component Semantic Versioning (`MAJOR.MINOR.PATCH`) identifiers.
For new versions, choose the appropriate change level based on compatibility
and functionality, and document the versioning convention consistently.
Do not rename historical tags.

A release normally needs only a short maintenance window of approximately
one or two days. Keep this focused on small fixes and cleanup rather than
introducing substantial new features.

### Release checklist

1. Review outstanding issues, resolve suitable small bugs, and record any
   known limitations that remain open.
2. Run the project's regression tests and check the GitHub Actions results.
   Run additional scientific or HPC validation when the changes warrant it.
3. Check documentation, citation metadata, and release-relevant version
   references. Confirm that the intended commit is ready for publication.
4. Select a version and create an annotated Git tag on the reviewed commit,
   for example `git tag -a vX.Y -m "Release vX.Y"`. Push the specific tag
   with `git push origin vX.Y`.
5. Summarize changes since the previous tag with
   `git log <previous-tag>..vX.Y --oneline`. Publish a GitHub Release
   from that tag with concise release notes and known limitations.
6. Check that Zenodo has archived the new release and provides a
   version-specific DOI. Use the project's concept DOI for references to
   the software across versions.
7. Register the new release in the JuSER publication database if applicable.

Release publication is a manual decision. A scheduled release window does
not require publishing code that is not ready.

## Cleaning the git repository

Always test the following procedures on a fresh copy of the
repository!

You can completely remove files from the Git repository using
`git-filter-repo`: https://github.com/newren/git-filter-repo

Analyze the current git repository:

```
    git filter-repo --analyze
    cd .git/filter-repo/analysis/
    ls
```

Remove all files except for those that currently exist:

```
    git ls-files >../paths-i-want.txt
    git filter-repo --paths-from-file ../paths-i-want.txt
```

## Further reading

- [Git Handbook](https://guides.github.com/introduction/git-handbook)

- [GitHub Git Cheat Sheet](https://training.github.com/downloads/github-git-cheat-sheet.pdf)
