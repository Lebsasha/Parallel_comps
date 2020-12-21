cmake_minimum_required(VERSION 3.16)
project(Lab_4)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS " -Weffc++ -Wextra -fopenmp ") #-Weffc++ -pedantic -Wextra

add_executable(Lab_4 main.cpp)