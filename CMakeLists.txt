SET(ARMADILLO_INCLUDE_DIRS include/armadillo-12.8.1/include/armadillo)

cmake_minimum_required(VERSION 3.23)
project(Proyecto)

set(CMAKE_CXX_STANDARD 23)

add_executable(Proyecto main.cpp src/tests.cpp src/Matrix.cpp include/Matrix.h src/EKF_Tests.cpp include/R_z.h src/R_z.cpp include/R_y.h src/R_y.cpp include/R_x.h src/R_x.cpp include/AccelPointMass.h src/AccelPointMass.cpp include/AzElPa.h src/AzElPa.cpp include/sign_.h src/sign_.cpp src/Cheb3D.cpp include/Cheb3D.h include/EccAnom.h src/EccAnom.cpp include/Frac.h src/Frac.cpp include/Geodetic.h src/SAT_const.cpp src/Geodetic.cpp include/SAT_const.h include/Global.h src/EKF_GEOS3.cpp src/Mjday.cpp include/Mjday.h src/Position.cpp include/Position.h src/Mjday_TDB.cpp include/Mjday_TDB.h src/timediff.cpp include/timediff.h src/unit.cpp include/unit.h src/NutAngles.cpp include/NutAngles.h src/MeanObliquity.cpp include/MeanObliquity.h src/Legendre.cpp include/Legendre.h src/IERS.cpp include/IERS.h src/TimeUpdate.cpp include/TimeUpdate.h src/NutMatrix.cpp include/NutMatrix.h src/PoleMatrix.cpp include/PoleMatrix.h src/PrecMatrix.cpp include/PrecMatrix.h src/angl.cpp include/angl.h src/dotProduct.cpp include/dotProduct.h src/MeasUpdate.cpp include/MeasUpdate.h src/doubler.cpp include/doubler.h src/crossProduct.cpp include/crossProduct.h src/gibbs.cpp include/gibbs.h src/gmst.cpp include/gmst.h src/elements.cpp include/elements.h src/EqnEquinox.cpp include/EqnEquinox.h src/G_AccelHarmonic.cpp include/G_AccelHarmonic.h src/gast.cpp include/gast.h src/GHAMatrix.cpp include/GHAMatrix.h src/hgibbs.cpp include/hgibbs.h src/LTC.cpp include/LTC.h src/VarEqn.cpp include/VarEqn.h src/AccelHarmonic.cpp include/AccelHarmonic.h src/Accel.cpp include/Accel.h src/JPL_Eph_DE430.cpp include/JPL_Eph_DE430.h src/anglesdr.cpp include/anglesdr.h src/anglesg.cpp include/anglesg.h src/largestRoot.cpp include/largestRoot.h include/types.h src/DEInteg.cpp include/DEInteg.h src/realmod.cpp include/realmod.h include/EFK_Tests.h)


