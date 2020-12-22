cmake_minimum_required(VERSION 3.16)
project(Lab_4)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS " -Weffc++ -Wextra -fopenmp ") #-Weffc++ -pedantic -Wextra

#//    180 200 220 240
#//    404 456 508 560
#//    628 712 796 880
#//    852 968 1084 1200


add_executable(Lab_4 main.cpp)