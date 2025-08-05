#!/bin/bash

set -e   # exit on any command failure

lode_png=src/lodepng.cpp
lode_png_out=src/lodepng.o
if [ ! -f "$lode_png_out" ]; then
   g++ -c "$lode_png" -o "$lode_png_out"
   if [ $? -ne 0 ]; then
      echo "Compiling "$lode_png" failed"
      exit 1
   fi
fi

utils_path=src/utils.cpp
utils_out=src/utils.o
g++ -c "$utils_path" -lfmt -o "$utils_out"
if [ $? -ne 0 ]; then
   echo "Compiling $utils_path failed"
   exit 2
fi

bidiagonalization_path=src/bidiagonalization.cpp
bidiagonalization_out=src/bidiagonalization.o
g++ "-I/usr/bin/eigen" -c "$bidiagonalization_path" -o "$bidiagonalization_out"
if [ $? -ne 0 ]; then
   echo "Compiling $bidiagonalization_path failed"
   exit 3 
fi

golub_kahan_path=src/golubKahan.cpp
golub_kahan_out=src/golubKahan.o
g++ "-I/usr/bin/eigen" -c "$golub_kahan_path" -o "$golub_kahan_out"
if [ $? -ne 0 ]; then
   echo "Compiling $golub_kahan_path failed"
   exit 4
fi

tests_path=src/tests.cpp
tests_out=src/tests.o
g++ "-I/usr/bin/eigen" -c "$tests_path" -o "$tests_out"
if [ $? -ne 0 ]; then
   echo "Compiling $tests_path failed"
   exit 5
fi

main=src/main.cpp
main_out=src/main.o
g++ "-I/usr/bin/eigen" -c $main -o "$main_out" 
if [ $? -ne 0 ]; then
   echo "Compiling $main failed"  
   exit 6
fi

g++ "$main_out" "$lode_png_out" "$utils_out" "$bidiagonalization_out" "$golub_kahan_out" "$tests_out" -lfmt -o main.exe
echo "main.cpp linked and compiled successfully"
./main.exe
