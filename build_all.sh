#!/bin/bash

for f in src/*.py; do
    echo "Building $f"

    docker run -v .:/host -w /host/src firedrakeproject/firedrake:latest python3 %%f
done