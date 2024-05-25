#include "../include/VarEqn.h"

extern Matrix eopdata;
extern aux AuxParam;

/**
 * @brief Computes the time derivative of the combined state vector and state transition matrix.
 *
 * This function calculates the time derivative of the combined state vector and state transition matrix
 * for the given time 'x' and the provided state vector 'yPhi'. It utilizes various auxiliary parameters
 * and functions to perform the computations.
 *
 * @param x Time at which to compute the time derivative, given in seconds.
 * @param yPhi Pointer to the Matrix object representing the combined state vector and state transition matrix.
 *             The state vector should be arranged as [r, v, Phi], where 'r' is the position vector,
 *             'v' is the velocity vector, and 'Phi' is the state transition matrix.
 * @return Matrix A Matrix object containing the time derivative of the combined state vector and state transition matrix.
 *
 * @note The function internally uses various auxiliary parameters and functions, such as IERS(), timediff(),
 *       PrecMatrix(), NutMatrix(), PoleMatrix(), GHAMatrix(), AccelHarmonic(), and G_AccelHarmonic(), to compute
 *       the required values and matrices.
 */
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

    Matrix yInput(yPhi->getColumnas(), yPhi->getFilas());
    yInput = yPhi->transpose();

    IERS(&eopdata, AuxParam.Mjd_UTC, 'l', x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);

    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;

    timediff(UT1_UTC,TAI_UTC, UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);

    double Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    // Transformation matrix
    Matrix P = PrecMatrix(MJD_J2000,AuxParam.Mjd_TT + x/86400.0);
    Matrix N = NutMatrix(AuxParam.Mjd_TT + x/86400.0);
    Matrix T = N * P;

    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    // State vector components
    Matrix r(3,1);
    r(1,1) = yInput(1,1);
    r(2,1) = yInput(2,1);
    r(3,1) = yInput(3,1);


    Matrix v(3,1);
    v(1,1) = yInput(4,1);
    v(2,1) = yInput(5,1);
    v(3,1) = yInput(6,1);


    Matrix Phi(6,6);



    // State transition matrix
    for(int j = 1; j <= 6; j++){
        for(int i = 1; i <= 6; i++){
            Phi(i, j) = yInput(6*j+i,1);
        }
    }

    // Acceleration and gradient
    Matrix a(3,1);
    a = AccelHarmonic ( &r, &E, AuxParam.n, AuxParam.m);
    Matrix G (3,3);
    G = G_AccelHarmonic ( &r, &E, AuxParam.n, AuxParam.m );

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

    Matrix Phip(6,6);
    Phip = dfdy*Phi;

    // Derivative of combined state vector and state transition matrix
    for(int i = 1; i <= 3; i++) {
        yPhip(i,1) = v(i,1);                 // dr/dt(i)
        yPhip(i + 3,1) = a(i,1);                 // dv/dt(i)
    }

    for(int i = 1; i <= 6; i++) {
        for (int j = 1; j <= 6; j++) {
            yPhip(6 * j + i, 1) = Phip(i, j);     //dPhi/dt(i,j)
        }
    }

    return yPhip.transpose();
}

Matrix VarEqnOUT(double x, Matrix& yPhi){
    return VarEqn(x, &yPhi);
}