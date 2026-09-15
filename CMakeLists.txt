cmake_minimum_required(VERSION 3.20)
project(navier_stokes_blowup_reproduction LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
enable_testing()
function(enable_project_warnings target)
 target_compile_options(${target} PRIVATE $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>)
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
add_executable(appendix_a3_diagnostics src/appendix_a3_diagnostics.cpp)
target_include_directories(appendix_a3_diagnostics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
enable_project_warnings(appendix_a3_diagnostics)
add_executable(appendix_a3_global_diagnostics src/appendix_a3_global_diagnostics.cpp)
target_include_directories(appendix_a3_global_diagnostics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
enable_project_warnings(appendix_a3_global_diagnostics)
add_executable(appendix_a4_pressure_diagnostics src/appendix_a4_pressure_diagnostics.cpp)
target_include_directories(appendix_a4_pressure_diagnostics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
enable_project_warnings(appendix_a4_pressure_diagnostics)
add_executable(appendix_a6_a7_diagnostics src/appendix_a6_a7_diagnostics.cpp)
target_include_directories(appendix_a6_a7_diagnostics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
enable_project_warnings(appendix_a6_a7_diagnostics)
add_executable(fftw_taylor_green src/fftw_taylor_green.cpp)
target_include_directories(fftw_taylor_green PRIVATE ${FFTW3_INCLUDE_DIR})
target_link_libraries(fftw_taylor_green PRIVATE ${FFTW3_LIBRARY})
enable_project_warnings(fftw_taylor_green)
add_executable(taylor_green_snapshots src/taylor_green_snapshots.cpp)
target_include_directories(taylor_green_snapshots PRIVATE ${FFTW3_INCLUDE_DIR})
target_link_libraries(taylor_green_snapshots PRIVATE ${FFTW3_LIBRARY})
enable_project_warnings(taylor_green_snapshots)
add_executable(fftw_taylor_green_3d src/fftw_taylor_green_3d.cpp)
target_include_directories(fftw_taylor_green_3d PRIVATE ${FFTW3_INCLUDE_DIR})
target_link_libraries(fftw_taylor_green_3d PRIVATE ${FFTW3_LIBRARY})
enable_project_warnings(fftw_taylor_green_3d)

set(NS_TESTS
 core_scaling:test_core_scaling.cpp
 similarity_coordinates:similarity_coordinates.cpp
 similarity_operators:similarity_operators.cpp
 leading_profile_relations:leading_profile_relations.cpp
 leading_stress:leading_stress.cpp
 leading_field_axis:leading_field_axis.cpp
 moment_corrections:moment_corrections.cpp
 profile_moment_patch:profile_moment_patch.cpp
 outer_profile_stages:outer_profile_stages.cpp
 appendix_a3_closure:appendix_a3_closure.cpp
 appendix_a3_global:appendix_a3_global.cpp
 appendix_a4_pressure:appendix_a4_pressure.cpp
 appendix_a6_heat:appendix_a6_heat.cpp
 appendix_a7_heat_compensation:appendix_a7_heat_compensation.cpp
 appendix_a7_full_schedule:appendix_a7_full_schedule.cpp
 appendix_a7_schedule_scales:appendix_a7_schedule_scales.cpp
 appendix_a7_log_discrepancy:appendix_a7_log_discrepancy.cpp
 appendix_a7_log_compensation:appendix_a7_log_compensation.cpp
 appendix_a7_full_recovery:appendix_a7_full_recovery.cpp
 appendix_b_axis:appendix_b_axis.cpp
 appendix_b_endpoint_interval:appendix_b_endpoint_interval.cpp
 appendix_b_Z_interval:appendix_b_Z_interval.cpp
 appendix_b_coefficient_space:appendix_b_coefficient_space.cpp
 appendix_b_contraction:appendix_b_contraction.cpp
 appendix_b_nonlinear_remainder:appendix_b_nonlinear_remainder.cpp
 appendix_b_proposition_b2_certificate:appendix_b_proposition_b2_certificate.cpp
 appendix_b_proposition_b3_certificate:appendix_b_proposition_b3_certificate.cpp
 divergence_free_projection:divergence_free_projection.cpp
 dealiasing_2_3:dealiasing_2_3.cpp
 viscous_decay:viscous_decay.cpp
 taylor_green:taylor_green_regression.cpp)
foreach(pair IN LISTS NS_TESTS)
 string(REPLACE ":" ";" fields "${pair}")
 list(GET fields 0 name)
 list(GET fields 1 file)
 add_executable(test_${name} tests/${file})
 if(NOT name STREQUAL "core_scaling" AND NOT name STREQUAL "divergence_free_projection" AND NOT name STREQUAL "dealiasing_2_3" AND NOT name STREQUAL "viscous_decay" AND NOT name STREQUAL "taylor_green")
  target_include_directories(test_${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
 endif()
 enable_project_warnings(test_${name})
 add_test(NAME ${name}_regression COMMAND test_${name})
endforeach()
add_test(NAME fourier_diff_1d_regression COMMAND fourier_diff_1d)
add_test(NAME fftw_taylor_green_smoke COMMAND fftw_taylor_green 6 0.02 0.01 0.1 ${CMAKE_CURRENT_BINARY_DIR}/fftw_taylor_green_smoke.csv)
add_test(NAME taylor_green_snapshot_smoke COMMAND taylor_green_snapshots 8 0.02 0.01 0.1 ${CMAKE_CURRENT_BINARY_DIR}/taylor_green_snapshots_smoke.csv)
add_test(NAME fftw_taylor_green_3d_smoke COMMAND fftw_taylor_green_3d 8 0.02 0.01 0.05 ${CMAKE_CURRENT_BINARY_DIR}/fftw_taylor_green_3d_smoke.csv 3)
