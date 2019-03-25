project(gpgpu)
cmake_minimum_required(VERSION 2.8)

# Compute capability. See -arch option in nvcc --help.
set (CC 35)

find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${CUDA_INCLUDE_DIRS}/../samples/common/inc)

set(CMAKE_CXX_FLAGS "-O3")
set(CUDA_PROPAGATE_HOST_FLAGS OFF)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-DSM_${CC};-arch=sm_${CC};-lineinfo;-Xcompiler=${CMAKE_CXX_FLAGS};-Xptxas=-v)

cuda_add_executable( outputCUDA wtime.h wtime.cu mmio.h mmio.cu print.h print.cu matrixPreprocessing.h matrixPreprocessing.cu matrixVector.h matrixVector.cu  test.h test.cu matrix_vector_CUDA.cu )
