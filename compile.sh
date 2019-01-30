#!/bin/bash

# if argument 1 is null string
if [ -z "$1" ] 
    then # set build type to Release
        echo "[WARNING] No argument supplied, expected CMAKE_BUILD_TYPE"
        echo "[WARNING] setting to Release"
        set -- "${@:1}" "Release"
        sleep 3
fi

echo "[BASH] rm -rf build/"
rm -rf build/
echo "[BASH] mkdir build"
mkdir build
echo "[BASH] cd build"
cd build
echo "[BASH] -DCMAKE_BUILD_TYPE=${1} -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++"
cmake .. -DCMAKE_BUILD_TYPE=$1 -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++
echo "[BASH] make -j"
make -j
# echo "[BASH] ctest"
# ctest
echo "[BASH] cd .."
cd ..
