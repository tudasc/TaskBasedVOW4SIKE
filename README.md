# vOW4SIKE: Task-based Implementation

Implementation accompanies the paper titled "Task-based Parallelization Approach for Attacking the Supersingular Isogeny Path Problem".
Our implementation is based on the implementation of Costello et al. which is published at [link](https://github.com/microsoft/vOW4SIKE).

**IMPORTANT:** 
All libraries (aes, sha3, xxhash, P128, vow_sike) which are found in folder `libs/` need to be compiled with the Intel Compiler in version as specified in the paper. 
For more details, see the attached script `script.sh`.
The `ml intel` command to load an Intel C compiler in Slurm can be commented out. 
Make sure that you use the Intel compiler and enable Intel TBB before compiling the code.

**Pre-defined macros:**
- {SERIAL_VER | SECOND_VER | TBB_VER}=ON/OFF
- BIG_MEM_FOR_OUTPUT=ON/OFF
- MITM=ON/OFF
- {p32 | p36 | p40 | p44 | p52 | p56}=ON/OFF
- PRECOMP_OUTPUT_TO_FILE=ON/OFF

## Build Instructions
### Step 1: Clone the xxHash library
In directory `libs` (in the main directory), run:
```
git clone https://github.com/Cyan4973/xxHash
```
to clone the xxHash library.

The hash function libraries xxhash, AES and SHA3 are free to use and they are originally used in the implementation of Costello et al.

Make sure that the `libs` directory contains the following subdirectories: `aes`, `P128`, `sha3`, `vow_sike`, `xxHash`.

### Step 2: Create build directories
Create a `build` directory inside each aforementioned subdirectory except for the `xxHash`.

### Step 3: Build libraries
- To build `aes`, `sha3` (and also for `vow_sike`, `P128`), run:
```
cmake ..
make
```
in the corresponding `build` directory. 

- To build `xxhash`, run:
`make` inside the `xxhash` directory.

- Assume that in the next steps we compile for the indexed-array approach and the prime characteristic `p32` and thus we specify the macros in CMake with:
```
cmake -DTBB_VER=ON -Dp32=ON -DCMAKE_BUILD_TYPE=Release ..
```
**Notes:** To switch to the Intel TBB approach, replace `TBB_VER` with `SECOND_VER`.  Similarly, to change the prime charateristic replace `p32` with one of the aforementioned primes.

To build `vow_sike` and `P128` run:
```
cmake -DTBB_VER=ON -Dp32=ON -DCMAKE_BUILD_TYPE=Release ..
make
```
**Notes:** Make sure that `vow_sike` is built before `P128`.

### Step 4: Build main project
Make a `build` directory inside the main directory. In `build` run:
```
cmake -DTBB_VER=ON -DBIG_MEM_FOR_OUTPUT=OFF -DMITM=OFF -Dp32=ON -DPRECOMP_OUTPUT_TO_FILE=OFF -DCMAKE_BUILD_TYPE=Release ..
make
```

## Execution Instructions
In general, run the code with:
`build/test_vOW_SIKE_128 -n num_threads -d delta`, where delta is the precomputation depth.

For more detailed usage and example scripts, refer to the provided `script.sh`.
