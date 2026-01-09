@echo off

docker run -v .:/host -w /host/src firedrakeproject/firedrake:latest python3 %*