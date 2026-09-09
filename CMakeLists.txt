cmake_minimum_required(VERSION 3.20)
project(navier_stokes_blowup_reproduction LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(core_scaling src/core_scaling.cpp)

target_compile_options(core_scaling PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
)
