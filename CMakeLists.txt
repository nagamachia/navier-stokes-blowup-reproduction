cmake_minimum_required(VERSION 3.20)
project(navier_stokes_blowup_reproduction LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

function(enable_project_warnings target)
    target_compile_options(${target} PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
    )
endfunction()

add_executable(core_scaling src/core_scaling.cpp)
enable_project_warnings(core_scaling)

add_executable(fourier_diff_1d src/fourier_diff_1d.cpp)
enable_project_warnings(fourier_diff_1d)

add_executable(test_core_scaling tests/test_core_scaling.cpp)
enable_project_warnings(test_core_scaling)

add_test(NAME core_scaling_regression COMMAND test_core_scaling)
add_test(NAME fourier_diff_1d_regression COMMAND fourier_diff_1d)
