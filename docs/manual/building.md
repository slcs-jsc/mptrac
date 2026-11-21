# Building from source

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

