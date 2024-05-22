#include "../include/Accel.h"

extern Matrix eopdata;
extern aux AuxParam;

//Y = (1,6)
//return (1,6)
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

    IERS(&eopdata, AuxParam.Mjd_UTC + x/86400,'l', x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);

    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;

    timediff(UT1_UTC,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);

    double Mjd_UT1 = AuxParam.Mjd_UTC + x/86400 + UT1_UTC/86400;
    double Mjd_TT = AuxParam.Mjd_UTC + x/86400 + TT_UTC/86400;

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
    Y1(1,1) = (*Y)(1,1);
    Y1(2,1) = (*Y)(1,2);
    Y1(3,1) = (*Y)(1,3);
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
    dY(1,1) = (*Y)(1,4);
    dY(1,2) = (*Y)(1,5);
    dY(1,3) = (*Y)(1,6);


    dY(1,4) = a(1,1);
    dY(1,5) = a(1,2);
    dY(1,6) = a(1,3);

    return dY;
}

//Y(6, 1)
//return (6,1)
Matrix AccelT(double x, Matrix* Y){
    Matrix YT(1,6);
    YT(1,1) = (*Y)(1,1);
    YT(1,2) = (*Y)(2,1);
    YT(1,3) = (*Y)(3,1);
    YT(1,4) = (*Y)(4,1);
    YT(1,5) = (*Y)(5,1);
    YT(1,6) = (*Y)(6,1);

    Matrix resT(1,6);
    resT = AccelT(x, &YT);

    Matrix res(6,1);
    res(1,1) = resT(1,1);
    res(2,1) = resT(1,2);
    res(3,1) = resT(1,3);
    res(4,1) = resT(1,4);
    res(5,1) = resT(1,5);
    res(6,1) = resT(1,6);

    return res;
}