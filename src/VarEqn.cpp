#include "../include/VarEqn.h"

Matrix VarEqn(double x, Matrix* yPhi){
    double x_pole;
    double y_pole;
    double UT1_UTC;
    double LOD;
    double dpsi;
    double deps;
    double dx_pole;
    double dy_pole;
    double TAI_UTC;

    IERS(Global::eopdata, Global::AuxParam.Mjd_UTC, 'l', x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);

    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;

    timediff(UT1_UTC,TAI_UTC, UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);

    double Mjd_UT1 = Global::AuxParam.Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    // Transformation matrix
    Matrix P = PrecMatrix(MJD_J2000,Global::AuxParam.Mjd_TT + x/86400);
    Matrix N = NutMatrix(Global::AuxParam.Mjd_TT + x/86400);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    // State vector components
    Matrix r(1,3);
    r(1,1) = (*yPhi)(1,1);
    r(1,2) = (*yPhi)(1,2);
    r(1,3) = (*yPhi)(1,3);


    Matrix v(1,3);
    v(1,1) = (*yPhi)(1,4);
    v(1,2) = (*yPhi)(1,5);
    v(1,3) = (*yPhi)(1,6);


    Matrix Phi(6,6);

    // State transition matrix
    for(int j = 1; j <= 6; j++){
        for(int i = 1; i <= 6; i++){
            Phi(i, j) = (*yPhi)(1,6*j+i);
        }
    }

    // Acceleration and gradient
    Matrix a = AccelHarmonic ( &r, &E, Global::AuxParam.n, Global::AuxParam.m);
    Matrix G = G_AccelHarmonic ( &r, &E, Global::AuxParam.n, Global::AuxParam.m );

    // Time derivative of state transition matrix
    Matrix yPhip(42,1);
    Matrix dfdy(6,6);

    for(int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            dfdy(i, j) = 0.0;                 // dv/dr(i,j)
            dfdy(i + 3, j) = G(i, j);            // da/dr(i,j)

            if (i == j) {
                dfdy(i, j + 3) = 1;
            } else {
                dfdy(i, j + 3) = 0;             // dv/dv(i,j)
            }

            dfdy(i + 3, j + 3) = 0.0;             // da/dv(i,j)
        }
    }

    Matrix Phip = dfdy*Phi;

    // Derivative of combined state vector and state transition matrix
    for(int i = 1; i <= 3; i++) {
        yPhip(i,1) = v(1,i);                 // dr/dt(i)
        yPhip(i + 3, 1) = a(1,i);                 // dv/dt(i)
    }

    for(int i = 1; i <= 6; i++) {
        for (int j = 1; j <= 6; j++) {
            yPhip(6 * j + i, 1) = Phip(i, j);     //dPhi/dt(i,j)
        }
    }

    return yPhip;
}