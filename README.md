# vOW4SIKE miniApp

Implementation accompanies the paper titled "Task-based Parallelization Approach for Attacking the Supersingular Isogeny Path Problem".
Our implementation is based on the implementation of Costello et al. which is published at [link](https://github.com/microsoft/vOW4SIKE).

**IMPORTANT:** 
All libraries (aes, sha3, xxhash, P128, vow_sike) which are found in folder `libs/` need to be compiled with the Intel Compiler in version as specified in the paper. 
For more details, see the attached script `script.sh`.
The `ml intel` command to load an Intel C compiler in Slurm can be commented out. 
Make sure that you use the Intel compiler and enable Intel TBB before compiling the code.

The pre-defined macros are:
- {SERIAL_VER|SECOND_VER|TBB_VER}=ON/OFF
- BIG_MEM_FOR_OUTPUT=ON/OFF
- MITM
- {p32 | p36 | p40 | p44 | p52 | p56}
- PRECOMP_OUTPUT_TO_FILE

Set `main_dir` to the main directory where this README resides.

In folder libs (in the main directory), run:
```
git clone https://github.com/Cyan4973/xxHash xxhash && cd xxhash && make && cd .. &&
mkdir -p P128/build aes/build sha3/build vow_sike/build &&
cd sha3/build && cmake .. && make && cd .. &&
cd aes/build && cmake .. && make
```
These commands above clone and build the hash function libraries xxHash, AES and SHA3, which are free to use. 
These libraries are originally used in the implementation of Costello et al.

To compile, run the following commands in the main directory:
```
mkdir build
cd build
cmake -D{SERIAL_VER|SECOND_VER|TBB_VER}=ON -DCMAKE_BUILD_TYPE={Release|RelWithDebInfo} ..
```

Run with:
`build/test_vOW_SIKE_128 -n num_threads -d delta`, where delta is the precomputation depth.
