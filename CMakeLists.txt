cmake_minimum_required(VERSION 3.16) 
# set the project name and version
project(test_vOW_SIKE_128 LANGUAGES CXX)

find_package(OpenMP REQUIRED)
find_package(TBB REQUIRED)

# specify the C++ standard
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_subdirectory(src)

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


if(PRECOMP_OUTPUT_TO_FILE)
add_definitions(-DPRECOMP_OUTPUT_TO_FILE)
endif()

if(PRINT_CURVE_INFO)
add_definitions(-DPRINT_CURVE_INFO)
endif()

if(BIG_MEM_FOR_OUTPUT)
add_definitions(-DBIG_MEM_FOR_OUTPUT)
endif()

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

# add the executable
file(GLOB SOURCES "${test_vOW_SIKE_128_SOURCE_DIR}/src/*.cpp")
file(GLOB HEADER "${test_vOW_SIKE_128_SOURCE_DIR}/include/*.h")
add_executable(test_vOW_SIKE_128 ${SOURCES} ${HEADER})

find_library(SIDH_LIB sidh HINTS ${CMAKE_CURRENT_SOURCE_DIR}/../P128/build)
find_library(VOW_SIKE_LIB vow_sike HINTS ${CMAKE_CURRENT_SOURCE_DIR}/libs/vow_sike/build)

# -l options
target_link_libraries(test_vOW_SIKE_128 PRIVATE 
        OpenMP::OpenMP_CXX
        ${VOW_SIKE_LIB}
        m
        TBB::tbb
)
# -I options
target_compile_options(test_vOW_SIKE_128 PRIVATE 
        ${OpenMP_CXX_FLAGS}
        -Wfatal-errors -march=native)
        
target_include_directories(test_vOW_SIKE_128 PRIVATE 
            ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/libs/vow_sike/include
            ${CMAKE_CURRENT_SOURCE_DIR}/libs/P128/include
            ${CMAKE_CURRENT_SOURCE_DIR}/libs/aes/include
            ${CMAKE_CURRENT_SOURCE_DIR}/libs/sha3/include
            ${CMAKE_CURRENT_SOURCE_DIR}/libs/xxhash
)
