#ifndef PROYECTO_IERS_H
#define PROYECTO_IERS_H


#include <iostream>
#include <vector>
#include "Matrix.h"
#include "SAT_const.h"

int findMatchRow(Matrix a, int b);

//x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC
void IERS(Matrix* eop, double Mjd_UTC, char interp, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC);

//x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC
void IERS(Matrix* eop, double Mjd_UTC, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC);


#endif //PROYECTO_IERS_H
