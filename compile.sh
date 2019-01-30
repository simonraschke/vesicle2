#!/bin/bash

# if argument 1 is null string
if [ -z "$1" ] 
    then # set build type to Release
        echo "[WARNING] No argument supplied, expected CMAKE_BUILD_TYPE"
        echo "[WARNING] setting to Release"
        set -- "${@:1}" "Release"
        sleep 3
fi

if [ -z "$2" ] 
    then # set build type to Release
        echo "[WARNING] No argument supplied, expected CMAKE_C_COMPILER"
        echo "[WARNING] setting to clang"
        set -- "${@:2}" "clang"
        sleep 3
fi

if [ -z "$3" ] 
    then # set build type to Release
        echo "[WARNING] No argument supplied, expected CMAKE_CXX_COMPILER"
        echo "[WARNING] setting to clang++"
        set -- "${@:3}" "clang++"
        sleep 3
fi

echo "[BASH] rm -rf build/"
rm -rf build/
echo "[BASH] mkdir build"
mkdir build
echo "[BASH] cd build"
cd build
echo "[BASH] -DCMAKE_BUILD_TYPE=${1} -DCMAKE_C_COMPILER=${2} -DCMAKE_CXX_COMPILER=${3}"
cmake .. -DCMAKE_BUILD_TYPE=$1 -DCMAKE_C_COMPILER=$2 -DCMAKE_CXX_COMPILER=$3
echo "[BASH] make -j"
make -j
# echo "[BASH] ctest"
# ctest
echo "[BASH] cd .."
cd ..
