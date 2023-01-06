cmake_minimum_required(VERSION 3.16)
project(sha3 LANGUAGES C)
find_package(OpenMP REQUIRED)

add_library(sha3 SHARED
    fips202.c
)

target_compile_options(sha3 PRIVATE -fopenmp -Wfatal-errors -mavx2 -maes -msse2)
target_include_directories(sha3 
PUBLIC include
)
                       

