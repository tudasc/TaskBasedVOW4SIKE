cmake_minimum_required(VERSION 3.16)
project(vow_sike LANGUAGES C)

find_package(OpenMP REQUIRED)

file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c")
file(GLOB HEADER "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")

add_library(vow_sike SHARED
        ${SOURCES}
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

if(PRECOMP_OUTPUT_TO_FILE)
add_definitions(-DPRECOMP_OUTPUT_TO_FILE)
endif()

if(PRINT_CURVE_INFO)
add_definitions(-DPRINT_CURVE_INFO)
endif()

if(BIG_MEM_FOR_OUTPUT)
add_definitions(-DBIG_MEM_FOR_OUTPUT)
endif()

IF (THIRD_VER)
add_definitions(-DTASK_BASED_PRECOMP_THIRD_VER) 
ENDIF()

IF (SECOND_VER)
add_definitions(-DTASK_BASED_PRECOMP_SECOND_VER) 
ENDIF()

if (TBB_VER)
add_definitions(-DTASK_BASED_PRECOMP_TBB_VER)
endif()

if(SERIAL_VER)
add_definitions(-DSERIAL_VER)
endif()

if(MITM)
add_definitions(-DMITM)
endif()

add_definitions(-DVOW_SIKE)
add_definitions(-DP128)
add_definitions(-D__LINUX__)
add_definitions(-D_AMD64_)
add_definitions(-D_FAST_)
add_definitions(-D_MULX_)
add_definitions(-D_ADX_)

find_library(SIDH_LIB sidh HINTS ${CMAKE_CURRENT_SOURCE_DIR}/../P128/build)
find_library(AES_LIB aes HINTS ${CMAKE_CURRENT_SOURCE_DIR}/../aes/build)
find_library(XXHASH_LIB xxhash HINTS ${CMAKE_CURRENT_SOURCE_DIR}/../xxhash)
find_library(SHA3_LIB sha3 HINTS ${CMAKE_CURRENT_SOURCE_DIR}/../sha3/build)

target_compile_options(vow_sike PUBLIC 
-fopenmp
-Wfatal-errors -mavx2 -maes -msse2)

target_link_libraries(vow_sike PRIVATE
${XXHASH_LIB} ${AES_LIB} ${SIDH_LIB} ${SHA3_LIB})

target_include_directories(vow_sike 
    PUBLIC include ${CMAKE_CURRENT_SOURCE_DIR}/../P128/include ${CMAKE_CURRENT_SOURCE_DIR}/../aes/include ${CMAKE_CURRENT_SOURCE_DIR}/../xxhash/
)
                       

