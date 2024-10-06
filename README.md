# vOW4SIKE: Task-based Implementation

Implementation accompanies the paper titled "Task-based Parallelization Approach for Attacking the Supersingular Isogeny Path Problem".
Our implementation is based on the implementation of Costello et al. which is published at [link](https://github.com/microsoft/vOW4SIKE).

**IMPORTANT:** 
All libraries (aes, sha3, xxhash, P128, vow_sike) which are found in folder `libs/` need to be compiled with the Intel Compiler in version as specified in the paper. 
For more details, see the attached script `script.sh`.
The `ml intel` command to load an Intel C compiler in Slurm can be commented out. 
Make sure that you use the Intel compiler and enable Intel TBB before compiling the code.

Pre-defined macros in CMake are:
- {SERIAL_VER|SECOND_VER|TBB_VER}=ON/OFF
- BIG_MEM_FOR_OUTPUT=ON/OFF
- MITM=ON/OFF
- {p32 | p36 | p40 | p44 | p52 | p56}=ON/OFF
- PRECOMP_OUTPUT_TO_FILE=ON/OFF

## Build Commands
In directory `libs` (in the main directory), run:
```
git clone https://github.com/Cyan4973/xxHash
```
to clone the xxHash library.

The hash function libraries xxhash, AES and SHA3 are free to use and they are originally used in the implementation of Costello et al.

Make sure that the `libs` directory contains the following subdirectories: `aes`, `P128`, `sha3`, `vow_sike`, `xxHash`.
Create a `build` directory inside each aforementioned subdirectory except for the `xxHash`.

Assume that in the next steps we compile for the indexed-array approach and the prime characteristic `p32` and thus we specify the macros every time we run CMake with:
```
cmake -DTBB_VER=ON -Dp32=ON -DCMAKE_BUILD_TYPE=Release ..
```
**Notes:** To switch to the Intel TBB approach, replace `TBB_VER` with `SECOND_VER`.  Similarly, to change the prime charateristic replace `p32` with one of the aforementioned primes.

To build `aes`, `sha3` (and also `P128` and `vow_sike`), run:
```
cmake -DTBB_VER=ON -Dp32=ON -DCMAKE_BUILD_TYPE=Release ..
make
```
in the corresponding `build` directory.

To build `xxhash`, run:
`make` inside the `xxhash` directory.

## Execution
Run the given script with:
`bash script.sh`
in the main directory, which by default builds the program for the prime `p32`.
The script then starts the program with 96 threads and sets the precomputation depth to 12. 

In general, run the code with:
`build/test_vOW_SIKE_128 -n num_threads -d delta`, where delta is the precomputation depth.
