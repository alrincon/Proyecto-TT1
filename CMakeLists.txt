SET(ARMADILLO_INCLUDE_DIRS include/armadillo-12.8.1/include/armadillo)

cmake_minimum_required(VERSION 3.23)
project(Proyecto)

set(CMAKE_CXX_STANDARD 23)

add_executable(Proyecto main.cpp src/tests.cpp src/Matrix.cpp include/Matrix.h src/EKF_Tests.cpp include/R_z.h src/R_z.cpp include/R_y.h src/R_y.cpp include/R_x.h src/R_x.cpp include/AzElPa.h src/AzElPa.cpp)