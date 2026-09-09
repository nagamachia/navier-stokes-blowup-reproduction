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

find_path(FFTW3_INCLUDE_DIR fftw3.h REQUIRED)
find_library(FFTW3_LIBRARY NAMES fftw3 REQUIRED)

add_executable(core_scaling src/core_scaling.cpp)
enable_project_warnings(core_scaling)

add_executable(fourier_diff_1d src/fourier_diff_1d.cpp)
enable_project_warnings(fourier_diff_1d)

add_executable(similarity_coordinates src/similarity_coordinates.cpp)
target_include_directories(similarity_coordinates PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
enable_project_warnings(similarity_coordinates)

add_executable(fftw_taylor_green src/fftw_taylor_green.cpp)
target_include_directories(fftw_taylor_green PRIVATE ${FFTW3_INCLUDE_DIR})
target_link_libraries(fftw_taylor_green PRIVATE ${FFTW3_LIBRARY})
enable_project_warnings(fftw_taylor_green)

add_executable(taylor_green_snapshots src/taylor_green_snapshots.cpp)
target_include_directories(taylor_green_snapshots PRIVATE ${FFTW3_INCLUDE_DIR})
target_link_libraries(taylor_green_snapshots PRIVATE ${FFTW3_LIBRARY})
enable_project_warnings(taylor_green_snapshots)

add_executable(test_core_scaling tests/test_core_scaling.cpp)
enable_project_warnings(test_core_scaling)

add_executable(test_similarity_coordinates tests/similarity_coordinates.cpp)
target_include_directories(test_similarity_coordinates PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
enable_project_warnings(test_similarity_coordinates)

add_executable(test_similarity_operators tests/similarity_operators.cpp)
target_include_directories(test_similarity_operators PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
enable_project_warnings(test_similarity_operators)

add_executable(test_leading_profile_relations tests/leading_profile_relations.cpp)
target_include_directories(test_leading_profile_relations PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
enable_project_warnings(test_leading_profile_relations)

add_executable(test_leading_stress tests/leading_stress.cpp)
target_include_directories(test_leading_stress PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
enable_project_warnings(test_leading_stress)

add_executable(test_leading_field_axis tests/leading_field_axis.cpp)
target_include_directories(test_leading_field_axis PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
enable_project_warnings(test_leading_field_axis)

add_executable(test_divergence_free_projection tests/divergence_free_projection.cpp)
enable_project_warnings(test_divergence_free_projection)

add_executable(test_dealiasing_2_3 tests/dealiasing_2_3.cpp)
enable_project_warnings(test_dealiasing_2_3)

add_executable(test_viscous_decay tests/viscous_decay.cpp)
enable_project_warnings(test_viscous_decay)

add_executable(test_taylor_green tests/taylor_green_regression.cpp)
enable_project_warnings(test_taylor_green)

add_test(NAME core_scaling_regression COMMAND test_core_scaling)
add_test(NAME similarity_coordinates_regression COMMAND test_similarity_coordinates)
add_test(NAME similarity_operators_regression COMMAND test_similarity_operators)
add_test(NAME leading_profile_relations_regression COMMAND test_leading_profile_relations)
add_test(NAME leading_stress_regression COMMAND test_leading_stress)
add_test(NAME leading_field_axis_regression COMMAND test_leading_field_axis)
add_test(NAME fourier_diff_1d_regression COMMAND fourier_diff_1d)
add_test(NAME divergence_free_projection_regression COMMAND test_divergence_free_projection)
add_test(NAME dealiasing_2_3_regression COMMAND test_dealiasing_2_3)
add_test(NAME viscous_decay_regression COMMAND test_viscous_decay)
add_test(NAME taylor_green_regression COMMAND test_taylor_green)
add_test(
    NAME fftw_taylor_green_smoke
    COMMAND fftw_taylor_green 6 0.02 0.01 0.1 ${CMAKE_CURRENT_BINARY_DIR}/fftw_taylor_green_smoke.csv
)
add_test(
    NAME taylor_green_snapshot_smoke
    COMMAND taylor_green_snapshots 8 0.02 0.01 0.1 ${CMAKE_CURRENT_BINARY_DIR}/taylor_green_snapshots_smoke.csv
)
