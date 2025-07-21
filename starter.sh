#!/bin/bash

set -e   # exit on any command failure

lode_png=/home/gabriel/Documents/HolyC/SVD_image_denoising/src/lodepng.cpp
lode_png_out=/home/gabriel/Documents/HolyC/SVD_image_denoising/src/lodepng.o
if [ ! -f "$lode_png_out" ]; then
   g++ -c "$lode_png" -o "$lode_png_out"
   if [ $? -ne 0 ]; then
      echo "Compiling "$lode_png" failed"
      exit 1
   fi
fi

image_utils_path=/home/gabriel/Documents/HolyC/SVD_image_denoising/src/imageUtils.cpp
image_utils_out=/home/gabriel/Documents/HolyC/SVD_image_denoising/src/imageUtils.o
g++ -c "$image_utils_path" -o "$image_utils_out"
if [ $? -ne 0 ]; then
   echo "Compiling $image_utils_path failed"
   exit 2
fi

bidiagonalization_path=/home/gabriel/Documents/HolyC/SVD_image_denoising/src/bidiagonalization.cpp
bidiagonalization_out=/home/gabriel/Documents/HolyC/SVD_image_denoising/src/bidiagonalization.o
g++ "-I/usr/bin/eigen" -c "$bidiagonalization_path" -o "$bidiagonalization_out"
if [ $? -ne 0 ]; then
   echo "Compiling $bidiagonalization_path failed"
   exit 3 
fi

golub_kahan_path=/home/gabriel/Documents/HolyC/SVD_image_denoising/src/golub_kahan.cpp
golub_kahan_out=/home/gabriel/Documents/HolyC/SVD_image_denoising/src/golub_kahan.hpp
g++ "-I/usr/bin/eigen" -c "$golub_kahan_path" -o "$golub_kahan_out"
if [ $? -ne 0 ]; then
   echo "Compiling $golub_kahan_path failed"
   exit 3 
fi

main=/home/gabriel/Documents/HolyC/SVD_image_denoising/src/main.cpp
main_out=/home/gabriel/Documents/HolyC/SVD_image_denoising/src/main.o
g++ "-I/usr/bin/eigen" -c $main -lfmt -o "$main_out" 
if [ $? -ne 0 ]; then
   echo "Compiling $main failed"  
   exit 4
fi

g++ "$main_out" "$lode_png_out" "$image_utils_out" "$bidiagonalization_out" "$golub_kahan_out" -lfmt -o main.exe
echo "main.cpp linked and compiled successfully"
./main.exe
