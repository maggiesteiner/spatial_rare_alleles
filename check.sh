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
if ! python bp_sims/source/simulations.py --json_out /dev/null -s 0.1 -L 50 --sampling_scheme "uniform" --mu 1e-08 --dens 1 --n_side 1 --time_limit 10 -r 0.1; then
    echo "FAIL: simulations.py"
    exit 1
fi
