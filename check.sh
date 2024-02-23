#!/usr/bin/env bash

echo Testing...
if ! ./bp_sims/source/test.py; then
    echo "FAIL: tests"
    exit 1
fi
