#!/bin/bash
# bash script will run as in a sandbox, so cd here won't change the pwd outside
# remember to `chmod +x compile.sh` for this file or use `sh compile.sh`
[[ -d build ]] || mkdir build
# rm build/CMakeCache.txt # if the directory is moved or reconstructed, reset the cmake cache before make
cd ./build
cmake ..
make
./FEM/main # be careful about the relative path, input/output file are relative to where you call the executable (rather than where the `main` exists)
