SET(ARMADILLO_INCLUDE_DIRS include/armadillo-12.8.1/include/armadillo)

cmake_minimum_required(VERSION 3.23)
project(Proyecto)

set(CMAKE_CXX_STANDARD 23)

add_executable(Proyecto main.cpp src/tests.cpp src/Matrix.cpp include/Matrix.h src/EKF_Tests.cpp include/R_z.h src/R_z.cpp include/R_y.h src/R_y.cpp include/R_x.h src/R_x.cpp include/AccelPointMass.h src/AccelPointMass.cpp include/AzElPa.h src/AzElPa.cpp include/sign_.h src/sign_.cpp src/Cheb3D.cpp include/Cheb3D.h include/EccAnom.h src/EccAnom.cpp include/Frac.h src/Frac.cpp include/Geodetic.h src/SAT_const.cpp src/Geodetic.cpp include/SAT_const.h)


