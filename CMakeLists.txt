cmake_minimum_required(VERSION 3.16)
project(mnogopotoki)
set(CMAKE_CXX_STANDARD 20)
add_executable(mnogopotoki main.cpp)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")
set(CMAKE_C_COMPILER "/usr/bin/mpicc")
set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")
include_directories(/usr/include/mpi/)
#"/usr/include/mpich/mpi.h"
