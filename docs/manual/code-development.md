# Code development

Please follow these guidelines when developing new code for MPTRAC.

## Compilation

To remove old binaries, backup files, etc., and to clean up the source code directory, use:

```
    make clean
```

For faster compilation, you can use parallel compilation with `make`:

```
    make -j
```

You can compile only a single binary, e.g.:

```
    make trac
```

You can change Makefile variables at the command line:

```
    make DEFINES="-DNP=100000" GPU=1 STATIC=0
```

## Testing

Always run all the tests to verify the revised code:

```
    make check
```

Please do not execute tests in parallel (`-j` option of `make`), you might miss seeing failed tests.

You can run only selected tests, e.g.:

```
    make trac_test
```

If a test fails, please carefully compare the test results with the reference data. The Linux tool `xxdiff` can be used to compare test data and reference data number by number. The graphing utility `gnuplot` can be used to visualize any differences. The reference data of a test should only be updated if the new test results are considered correct.

On Linux, use the `indent` tool to apply a selected formatting style to the source code:

```
    make indent
```

Only use `indent` if the code has been compiled correctly. You may need to re-run the `indent` command 2-3 times to get proper results.

You can use static code analysis to automatically detect potential problems in the code:

```
    make cppcheck
    make lizard
```

You can perform a coverage analysis to determine which parts of the code are covered by tests:

```
    make COV=1
    make coverage
```

If you find that parts of the code are not covered by the coverage analysis, please consider adding a new test.

After committing revised code to the GitHub repository, please check the [GitHub Actions page](https://github.com/slcs-jsc/mptrac/actions) to see if the automated tests were successfully passed.

Please also check the [Codacy](https://app.codacy.com/gh/slcs-jsc/mptrac?utm_source=github.com&utm_medium=referral&utm_content=slcs-jsc/mptrac&utm_campaign=Badge_Grade_Settings) and [Codecov](https://codecov.io/gh/slcs-jsc/mptrac) websites for test results.

Please check the [nightly build website](https://datapub.fz-juelich.de/slcs/mptrac/nightly_builds) to see if the tests passed on the Juelich supercomputers, especially for the GPU code. User at the Juelich Supercomputing Centre can also check the test results at the [ESM Buildbot](https://esm-buildbot.fz-juelich.de).

## Documentation

To update the [Doxygen documentation](https://slcs-jsc.github.io/mptrac/) (HTML and pdf files) from the source code, use:

```
    make doc
```

Please update the [GitHub wiki pages](https://github.com/slcs-jsc/mptrac/wiki/) to describe all changes and new code.

Please update the [README file](https://github.com/slcs-jsc/mptrac/blob/master/README.md) as needed.

## Installation

To copy the executables from the source directory to the `DESTDIR` directory (default `../bin/`):

```
    make install
```

To remove the executables from the installation directory:

```
    make uninstall
```

To create a zip file of the current state of the repository, including the compiled binaries, run

```
    make dist
```

## Releases

A new release of MPTRAC is usually made every six months.

To create a new release, first define a version tag ("vX.Y")  in the local repository:

```
    gitk
```

Next, push the local tags to the remote repository on GitHub:

```
    git push --tags
```

Get a list of short log messages from the previous to the current version of the code:

```
    git log v1.1..v1.2 --oneline
```

Using the log messages, [draft a new release](https://github.com/slcs-jsc/mptrac/releases/new) on GitHub using the newly created tag.

Check the [Zenodo the website](https://doi.org/10.5281/zenodo.4400597) to publish the new release and to get a DOI.

Use the DOI to submit entry at [JuSER publication database](https://juser.fz-juelich.de/).

## Cleaning the git repository

Always test the following procedures on a fresh copy of the repository!

Files can be completely removed from the git repository using `git-filter-repo`: https://github.com/newren/git-filter-repo

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
