SET(ARMADILLO_INCLUDE_DIRS include/armadillo-12.8.1/include/armadillo)

cmake_minimum_required(VERSION 3.23)
project(Proyecto)

set(CMAKE_CXX_STANDARD 23)

add_executable(Proyecto main.cpp include/matriz.cpp include/matriz.h)