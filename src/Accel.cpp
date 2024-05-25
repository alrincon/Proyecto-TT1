#include "../include/Accel.h"

extern Matrix eopdata;
extern aux AuxParam;


//------------------------------------------------------------------------------
// Matrix Accel(double x, Matrix* Y)
//------------------------------------------------------------------------------
/**
 * Calculates the acceleration of a spacecraft based on various perturbations.
 * This function computes the acceleration due to harmonic gravity field, luni-solar
 * perturbations, and planetary perturbations.
 *
 * @param <x> time since epoch in seconds
 * @param <Y> pointer to a Matrix containing state vector [position; velocity]
 * @return a Matrix containing the derivative of the state vector [velocity; acceleration]
 * @exception none
 * @note caller is responsible for managing the memory of the input Matrix Y
 */
//------------------------------------------------------------------------------
Matrix Accel(double x, Matrix* Y){
    double x_pole;
    double y_pole;
    double UT1_UTC;
    double LOD;
    double dpsi;
    double deps;
    double dx_pole;
    double dy_pole;
    double TAI_UTC;

    Matrix yInput(Y->getColumnas(), Y->getFilas());
    yInput = Y->transpose();

    IERS(&eopdata, AuxParam.Mjd_UTC + x/86400.0,'l', x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);

    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;

    timediff(UT1_UTC,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);

    double Mjd_UT1 = AuxParam.Mjd_UTC + x/86400.0 + UT1_UTC/86400.0;
    double Mjd_TT = AuxParam.Mjd_UTC + x/86400.0 + TT_UTC/86400.0;

    Matrix P = PrecMatrix(MJD_J2000,Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    double MJD_TDB = Mjday_TDB(Mjd_TT);

    Matrix r_Mercury(1,3);
    Matrix r_Venus(1,3);
    Matrix r_Earth(1,3);
    Matrix r_Mars(1,3);
    Matrix r_Jupiter(1,3);
    Matrix r_Saturn(1,3);
    Matrix r_Uranus(1,3);
    Matrix r_Neptune(1,3);
    Matrix r_Pluto(1,3);
    Matrix r_Moon(1,3);
    Matrix r_Sun(1,3);

    JPL_Eph_DE430(MJD_TDB, r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, r_Neptune,r_Pluto,r_Moon,r_Sun);

    // Acceleration due to harmonic gravity field

    Matrix Y1(3,1);
    Y1(1,1) = yInput(1,1);
    Y1(2,1) = yInput(2,1);
    Y1(3,1) = yInput(3,1);

    Matrix a(3,1);
    a = AccelHarmonic(&Y1, &E, AuxParam.n, AuxParam.m);

    // Luni-solar perturbations
    if (AuxParam.sun > 0){
        a = a + AccelPointMassT(Y1,r_Sun,GM_Sun);
    }

    if (AuxParam.moon > 0) {
        a = a + AccelPointMassT(Y1,r_Moon,GM_Moon);
    }

    // Planetary perturbations
    if (AuxParam.planets > 0) {
        a = a + AccelPointMassT(Y1,r_Mercury,GM_Mercury);
        a = a + AccelPointMassT(Y1,r_Venus,GM_Venus);
        a = a + AccelPointMassT(Y1,r_Mars,GM_Mars);
        a = a + AccelPointMassT(Y1,r_Jupiter,GM_Jupiter);
        a = a + AccelPointMassT(Y1,r_Saturn,GM_Saturn);
        a = a + AccelPointMassT(Y1,r_Uranus,GM_Uranus);
        a = a + AccelPointMassT(Y1,r_Neptune,GM_Neptune);
        a = a + AccelPointMassT(Y1,r_Pluto,GM_Pluto);
    }

    Matrix dY(1,6);
    dY(1,1) = yInput(4,1);
    dY(1,2) = yInput(5,1);
    dY(1,3) = yInput(6,1);

    dY(1,4) = a(1,1);
    dY(1,5) = a(2,1);
    dY(1,6) = a(3,1);

    return dY;
}


//------------------------------------------------------------------------------
// Matrix AccelOUT(double x, Matrix& Y)
//------------------------------------------------------------------------------
/**
 * Wrapper function for Accel.
 * This function provides a more convenient interface for calling Accel with a Matrix reference.
 *
 * @param <x> time since epoch in seconds
 * @param <Y> reference to a Matrix containing state vector [position; velocity]
 * @return a Matrix containing the derivative of the state vector [velocity; acceleration]
 * @exception none
 * @note no special notes
 */
//------------------------------------------------------------------------------
Matrix AccelOUT(double x, Matrix& Y){
    return Accel(x, &Y);
}