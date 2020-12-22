cmake_minimum_required(VERSION 3.16)
project(Temp_mpi)

set(CMAKE_CXX_STANDARD 14)

include_directories(. /usr/include/x86_64-linux-gnu/mpich/)
set(CMAKE_C_COMPILER "/usr/bin/mpicc")
set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")
set(CMAKE_CXX_FLAGS " -Weffc++ -pedantic -Wextra ")
add_executable(Temp_mpi
        main.cpp)
