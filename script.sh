#!/bin/sh
#ml intel
#ml cmake
spack load intel-tbb
eval `spack load --sh   intel-tbb`

main_dir="./" 
output_dir="./"

#modulus_arr=(p32 p36 p40 p44 p48 p52 p56)
#d_arr=(14 16 18 20 22 24 26) # mitm
#d_arr=(12  14  16  18  20  22  24) #vow

modulus=p32
d_arr=(12)
repeat=1

for (( i=0; i < ${#d_arr[@]}; ++i ))   
do
d=${d_arr[$i]}

build_options="-DTBB_VER=ON -DBIG_MEM_FOR_OUTPUT=OFF -DMITM=OFF -D$modulus=ON -DPRECOMP_OUTPUT_TO_FILE=OFF -DCMAKE_BUILD_TYPE=Release"
cd $main_dir
cd libs
cd P128/build && rm -rf * && cmake $build_options .. && make
cd $main_dir
cd libs
cd vow_sike/build && rm -rf * && cmake $build_options .. && make
cd $main_dir
cd build && rm -rf * && cmake $build_options .. && make

for (( j=0; j < ${repeat}; ++j ))
do
./test_vOW_SIKE_128 -n 96 -d $d 
done
done
