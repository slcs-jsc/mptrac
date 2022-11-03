# Code development

Please follow these guidelines when developing new code for MPTRAC.

## Compilation

To remove old binaries, backup files, etc., and to clean the source code directory, use:

```
    make clean
```

For faster compilation, you can apply parallel compilation with `make`, e.g.:

```
    make -j4
```

You can compile only a single binary, e.g.:

```
    make trac
```

You can modify the variables of the Makefile at the command line, e.g.:

```
    make DEFINES="-DNP=100000" GPU=1 STATIC=0
```

## Testing

Always run all the tests to check the revised code:

```
    make check
```

Please do not execute tests in parallel (option `-j` of `make`), you may overlook failed tests in this case.

You can run only selected tests, e.g.:

```
    make trac_test
```

If a test fails, please carefully compare the test results with the reference data. The Linux tool `xxdiff` might be used to compare test data and reference data number by number. The graphing utility `gnuplot` may be used to visualize any differences. The reference data of a test should only be updated if the new test results are considered correct.

Use the `indent` tool on Linux to apply the given format style to the source code:

```
    make indent
```

Only use `indent`, when the code was compiled correctly. You may have to rerun the `indent` command 2-3 times to get proper results.

You can apply static code analysis to automatically identify potential issues with the code:

```
    make cppcheck
    make lizard
```

You can run a coverage analysis to identify the parts of code which are covered by tests:

```
    make COV=1
    make coverage
```

If you find parts of the code are not covered in the coverage analysis, please consider adding a new test.

After pushing revised code to the GitHub repository, please check the [GitHub Actions page](https://github.com/slcs-jsc/mptrac/actions) to see whether the automatic tests were successfully passed.

Please also check the [Codacy](https://app.codacy.com/gh/slcs-jsc/mptrac?utm_source=github.com&utm_medium=referral&utm_content=slcs-jsc/mptrac&utm_campaign=Badge_Grade_Settings) and [Codecov](https://codecov.io/gh/slcs-jsc/mptrac) websites to see the test results.

Please check the [nightly build website](https://datapub.fz-juelich.de/slcs/mptrac/nightly_builds) to see whether checks on the Juelich supercomputers were passed.

## Documentation

To update the [Doxygen documentation](https://slcs-jsc.github.io/mptrac/) (HTML and pdf files) from the source code:

```
    make doc
```

Please update the [GitHub wiki pages](https://github.com/slcs-jsc/mptrac/wiki/) to describe any changes and new code.

Please update the [README file](https://github.com/slcs-jsc/mptrac/blob/master/README.md) as needed.

## Installation

To copy the executables from the source directory to the directory `DESTDIR` (default `../bin/`):

```
    make install
```

To remove the executables from the installation directory:

```
    make uninstall
```

To create a zip file of the current state of the repository, including the compiled binaries:

```
    make dist
```

## Releases

A new release of MPTRAC is typically created every six months.

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

Check the [Zenodo the website](https://doi.org/10.5281/zenodo.4400597) to publish the new release and to receive a DOI.

## Cleaning the git repository

Always test the procedures described below on a fresh copy of the repository!

Files can be removed completely from the git repository by means of `git-filter-repo`: https://github.com/newren/git-filter-repo

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

[Git Handbook](https://guides.github.com/introduction/git-handbook)

[GitHub Git Cheat Sheet](https://training.github.com/downloads/github-git-cheat-sheet.pdf)
