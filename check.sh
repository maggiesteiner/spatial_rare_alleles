#!/usr/bin/env bash

echo Testing...
if ! ./bp_sims/source/test.py; then
    echo "FAIL: tests"
    exit 1
fi

echo Running mypy to check types...
if ! mypy --pretty bp_sims/source/; then
    echo "FAIL: types"
    exit 1
fi

echo Fuzzing simulations.py...
if ! python bp_sims/source/simulations.py --sampled_p_out /dev/null --zero_out /dev/null; then
    echo "FAIL: simulations.py"
    exit 1
fi
