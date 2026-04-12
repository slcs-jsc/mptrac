#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C
set -euo pipefail

# Setup...
trac=../../src

# Create directory...
rm -rf data && mkdir -p data

# Test common CLI behavior for every executable listed in src/Makefile.
# Keeping the binary list in the Makefile avoids a second list here that
# would need manual updates whenever a tool is added or removed.
error=0
while IFS= read -r tool ; do
    [ -n "$tool" ] || continue
    echo "checking $tool..."

    # Calling a tool without arguments must fail and print the standard
    # "missing arguments" diagnostic used throughout the command-line tools.
    noarg_out="data/${tool}_noarg.txt"
    if "$trac/$tool" > "$noarg_out" 2>&1 ; then
        echo "Missing-argument command succeeded: $tool"
        error=1
    else
        if ! grep -q "Missing or invalid command-line arguments." "$noarg_out"; then
            echo "Missing argument error text: $tool"
            error=1
        fi
    fi

    # Both documented help flags must succeed and print a usage section.
    for flag in -h --help ; do
        echo "  help flag $flag"
        outfile="data/${tool}_${flag#--}.txt"
        if ! "$trac/$tool" "$flag" > "$outfile" ; then
            echo "Help command failed: $tool $flag"
            error=1
            continue
        fi
        if ! grep -q "Usage:" "$outfile" ; then
            echo "Missing usage text: $tool $flag"
            error=1
        fi

        # Help should be discoverable even when users add extra arguments.
        # This keeps help behavior forgiving for quick terminal exploration.
        extra_out="data/${tool}_${flag#--}_extra.txt"
        if ! "$trac/$tool" "$flag" extra-arg > "$extra_out" ; then
            echo "Help command with extra arguments failed: $tool $flag"
            error=1
            continue
        fi
        if ! grep -q "Usage:" "$extra_out" ; then
            echo "Missing usage text with extra arguments: $tool $flag"
            error=1
        fi
    done

# Read the executable list from the Makefile so this test follows the same
# source of truth as the build and installation rules.
done < <(make --no-print-directory -C ../../src print-exc)

exit $error
