cmake_minimum_required(VERSION 3.16)
project(sidh LANGUAGES C ASM)

add_library(sidh STATIC
        src/AMD64/fp_x64_asm.S
        src/AMD64/fp_x64.c
        src/fpx.c
        src/ec_isogeny.c
        src/P128.c
)

if(p32)
add_definitions(-Dp_32_20)
endif()

if(p36)
add_definitions(-Dp_36_22)
endif()

if(p40)
add_definitions(-Dp_40_25)
endif()

if(p44)
add_definitions(-Dp_44_27)
endif()

if(p48)
add_definitions(-Dp_48_30)
endif()

if(p52)
add_definitions(-Dp_52_33)
endif()

if(p56)
add_definitions(-Dp_56_35)
endif()

add_definitions(-DP128)
add_definitions(-D__LINUX__)
add_definitions(-D_AMD64_)
add_definitions(-D_FAST_)
add_definitions(-D_MULX_)
add_definitions(-D_ADX_)

target_compile_options(sidh PRIVATE -fopenmp -fPIC -Wfatal-errors -mavx2 -maes -msse2)
target_include_directories(sidh 
    PUBLIC include
)

                       

