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
