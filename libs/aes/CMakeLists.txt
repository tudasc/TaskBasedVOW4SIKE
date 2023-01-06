cmake_minimum_required(VERSION 3.16)
project(aes LANGUAGES C)
find_package(OpenMP REQUIRED)

add_library(aes SHARED
    aes_c.c
    aes_ni.c
    aes.c
)

target_compile_options(aes PRIVATE -fopenmp -Wfatal-errors -mavx2 -maes -msse2)
target_include_directories(aes 
    PUBLIC include
    PRIVATE .                 
)
                       

