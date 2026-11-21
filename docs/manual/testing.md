# Testing and code quality

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

