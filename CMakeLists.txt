cmake_minimum_required(VERSION 3.13.4)

project(Capture LANGUAGES C CUDA)

# Build Configuration
find_package(CUDAToolkit REQUIRED)

set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CUDA_RUNTIME_LIBRARY Static)

# Set arch
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES "50;75;80")
endif()

if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 14)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

# Chris put all the headers in one place, so
include_directories(".")

# PSRDada
find_package(PSRDada)

# Thrust
find_package(Thrust REQUIRED CONFIG)
thrust_create_target(Thrust)

# CUDA Flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=-Wall ")

# All the executables
add_executable(dsaX_correlator_udpdb_thread dsaX_correlator_udpdb_thread.c)
add_executable(dsaX_spectrometer_udpdb_thread dsaX_spectrometer_udpdb_thread.c)
add_executable(dsaX_spectrometer_reorder dsaX_spectrometer_reorder.cu)
add_executable(dsaX_spectrometer dsaX_spectrometer.cu)
add_executable(dsaX_correlator_trigger dsaX_correlator_trigger.c)
add_executable(dsaX_dbdisk dsaX_dbdisk.c)
add_executable(dsaX_dbdisk_fb dsaX_dbdisk_fb.c)

# Includes
target_link_libraries(dsaX_correlator_udpdb_thread psrdada)
target_link_libraries(dsaX_spectrometer_udpdb_thread psrdada)
target_link_libraries(dsaX_spectrometer_reorder Thrust psrdada)
target_link_libraries(dsaX_spectrometer Thrust psrdada m)
target_link_libraries(dsaX_correlator_trigger psrdada m)
target_link_libraries(dsaX_dbdisk psrdada m)
target_link_libraries(dsaX_dbdisk_fb psrdada m)
