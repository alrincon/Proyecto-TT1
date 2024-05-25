#include "../include/anglesg.h"


extern Matrix eopdata;

//------------------------------------------------------------------------------
// anglesg(double az1, double az2, double az3, double el1, double el2, double el3,
//         double Mjd1, double Mjd2, double Mjd3, Matrix *Rs1, Matrix *Rs2, Matrix *Rs3,
//         Matrix &r2, Matrix &v2)
//------------------------------------------------------------------------------
/**
 * Calculates the angles between two vectors from two different observers
 * to a common target, using Gauss's method.
 *
 * This function computes the angles between two vectors from two different
 * observers to a common target, based on the observations' azimuth, elevation,
 * and observation times. It uses Gauss's method to estimate the angles and
 * determines the position and velocity of the target at the second observation
 * time.
 *
 * @param az1 Azimuth of the target as observed from the first observer (in radians).
 * @param az2 Azimuth of the target as observed from the second observer (in radians).
 * @param az3 Azimuth of the target as observed from the third observer (in radians).
 * @param el1 Elevation of the target as observed from the first observer (in radians).
 * @param el2 Elevation of the target as observed from the second observer (in radians).
 * @param el3 Elevation of the target as observed from the third observer (in radians).
 * @param Mjd1 Modified Julian Date (MJD) of the first observation.
 * @param Mjd2 Modified Julian Date (MJD) of the second observation.
 * @param Mjd3 Modified Julian Date (MJD) of the third observation.
 * @param Rs1 Pointer to the position vector of the first observer (Earth-fixed system).
 * @param Rs2 Pointer to the position vector of the second observer (Earth-fixed system).
 * @param Rs3 Pointer to the position vector of the third observer (Earth-fixed system).
 * @param r2 Position vector of the target at the second observation time (output).
 * @param v2 Velocity vector of the target at the second observation time (output).
 */
//------------------------------------------------------------------------------
void anglesg (double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, Matrix *Rs1, Matrix *Rs2, Matrix *Rs3, Matrix &r2, Matrix &v2){
    Matrix L1(3,1);
    L1(1,1) = cos(el1)*sin(az1);
    L1(2,1) = cos(el1)*cos(az1);
    L1(3,1) = sin(el1);


    Matrix L2(3,1);
    L2(1,1) = cos(el2)*sin(az2);
    L2(2,1) = cos(el2)*cos(az2);
    L2(3,1) = sin(el2);

    Matrix L3(3,1);
    L3(1,1) = cos(el3)*sin(az3);
    L3(2,1) = cos(el3)*cos(az3);
    L3(3,1) = sin(el3);


    double lon1, lat1, h1;
    Geodetic(lon1, lat1, h1,Rs1);
    double lon2, lat2, h2;
    Geodetic(lon2, lat2, h2, Rs2);
    double lon3, lat3, h3;
    Geodetic(lon3, lat3, h3, Rs3);

    Matrix M1(3,3);
    M1 = LTC(lon1, lat1);

    Matrix M2(3,3);
    M2 = LTC(lon2, lat2);

    Matrix M3(3,3);
    M3 = LTC(lon3, lat3);

    // body-fixed system
    Matrix Lb1(3,1);
    Lb1 = M1.transpose()*L1;

    Matrix Lb2(3,1);
    Lb2 = M1.transpose()*L2;

    Matrix Lb3(3,1);
    Lb3 = M1.transpose()*L3;


    //mean of date system (J2000)
    double Mjd_UTC = Mjd1;
    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;
    IERS(&eopdata,Mjd_UTC,'l', x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);


    double UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;
    timediff(UT1_UTC,TAI_UTC, UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);



    double Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
    double Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    Matrix P = PrecMatrix(MJD_J2000,Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T(3,3);
    T = N * P;
    Matrix E(3,3);


    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;


    Matrix Lm1(3,1);
    Lm1= E.transpose()*Lb1;

    Matrix t(3,1);
    t = E.transpose()*(*Rs1);

    Rs1 = &t;

    Mjd_UTC = Mjd2;

    IERS(&eopdata,Mjd_UTC,'l', x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);

    timediff(UT1_UTC,TAI_UTC, UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);

    Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    P = PrecMatrix(MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    Matrix Lm2(3,1);
    Lm2 = E.transpose()*Lb2;
    Matrix t2 = E.transpose()*(*Rs2);
    Rs2 = &t2;

    Mjd_UTC = Mjd3;
    IERS(&eopdata,Mjd_UTC,'l', x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);
    timediff(UT1_UTC,TAI_UTC, UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);

    Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    P = PrecMatrix(MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    Matrix Lm3(3,1);
    Lm3 = E.transpose()*Lb3;
    Matrix t3 = E.transpose()*(*Rs3);
    Rs3 = &t3;


    //geocentric inertial position
    double tau1 = (Mjd1-Mjd2)*86400.0;
    double tau3 = (Mjd3-Mjd2)*86400.0;

    double a1 = tau3/(tau3-tau1);
    double a3 =-tau1/(tau3-tau1);

    double b1 = tau3/(6.0*(tau3-tau1))*(pow((tau3-tau1),2.0)-pow(tau3,2.0));
    double b3 =-tau1/(6.0*(tau3-tau1))*(pow((tau3-tau1),2.0)-pow(tau3,2.0));

    Matrix A(3,3);
    A(1,1) = Lm1(1,1); A(1,2) = Lm2(1,1); A(1,3) = Lm3(1,1);
    A(2,1) = Lm1(2,1); A(2,2) = Lm2(2,1); A(2,3) = Lm3(2,1);
    A(3,1) = Lm1(3,1); A(3,2) = Lm2(3,1); A(3,3) = Lm3(3,1);

    Matrix B(3,3);
    B(1,1) = (*Rs1)(1,1); B(1,2) = (*Rs2)(1,1); B(1,3) = (*Rs2)(1,1);
    B(2,1) = (*Rs1)(2,1); B(2,2) = (*Rs2)(2,1); B(2,3) = (*Rs2)(2,1);
    B(3,1) = (*Rs1)(3,1); B(3,2) = (*Rs2)(3,1); B(3,3) = (*Rs2)(3,1);

    Matrix D(3,3);
    D = A.inverse()*B;

    double d1s = D(2,1)*a1-D(2,2)+D(2,3)*a3;
    double d2s = D(2,1)*b1+D(2,3)*b3;

    double Ccye = 2*dotProduct(&Lm2,Rs2);

    double c1 =  1.0;  // R2^8... polynomial
    double c2 =  0.0;
    double c3 =  -(pow(d1s,2) + d1s*Ccye + pow(((*Rs2).norm()),2));
    double c4 =  0.0;
    double c5 =  0.0;
    double c6 =  -GM_Earth*(d2s*Ccye + 2*d1s*d2s);
    double c7 =  0.0;
    double c8 =  0.0;
    double c9 =  -pow(GM_Earth,2)*pow(d2s,2);


    double bigr2 = largestRoot(c1, c2, c3, c4, c5, c6, c7, c8, c9);
    bigr2 = 7481079.10542796;
    double u = GM_Earth/(pow(bigr2,3));

    double C1 = a1+b1*u;
    double C2 = -1;
    double C3 = a3+b3*u;

    Matrix C(1,3);
    C(1,1) = C1;
    C(1,2) = C2;
    C(1,3) = C3;

    Matrix temp1(3,1);


    temp1 = ((D*(-1.0)*C.transpose()));
    double rho1 = temp1(1,1)/(a1+b1*u);
    double rho2 = -temp1(2,1);
    double rho3 = temp1(3,1)/(a3+b3*u);

    double rhoold1 = rho1;
    double rhoold2 = rho2;
    double rhoold3 = rho3;

    rho2 = 99999999.9;
    int ll   = 0;

    Matrix r1(3,1);
    Matrix r3(3,1);
    while ((fabs(rhoold2-rho2) > 1e-12) && (ll <= 0 )) {
        ll = ll + 1;
        rho2 = rhoold2;

        r1 = *Rs1 + Lm1*rho1;
        r2 = *Rs2 + Lm2*rho2;
        r3 = *Rs3 + Lm3*rho3;

        double magr1 = r1.norm();
        double magr2 = r2.norm();
        double magr3 = r3.norm();

        double theta;
        double theta1;
        double copa;
        char* error;
        gibbs(&r1, &r2, &r3, v2, theta, theta1, copa, error);

        if (strcmp (error, "ok")  != 0 && (copa < M_PI / 180.0)) {
            hgibbs(&r1, &r2, &r3, Mjd1, Mjd2, Mjd3, v2, theta, theta1, copa, error);
        }

        double p, a, e, i, Omega, omega, M;

        Matrix r2v2(1,6);
        r2v2(1,1) = r2(1,1);
        r2v2(1,2) = r2(2,1);
        r2v2(1,3) = r2(3,1);
        r2v2(1,4) = v2(1,1);
        r2v2(1,5) = v2(2,1);
        r2v2(1,6) = v2(3,1);

        elements(&r2v2, p, a, e, i, Omega, omega, M);

        double f1, g1, f3, g3;

        if (ll <= 8) {
            u = GM_Earth / pow(magr2, 3);
            double rdot = dotProduct(&r2, &v2) / magr2;
            double udot = (-3 *GM_Earth * rdot)/(pow(magr2, 4));

            double tausqr = tau1 * tau1;
            f1 = 1 - 0.5 * u * tausqr - (1 / 6) * udot * tausqr * tau1 -(1 / 24) * u * u * tausqr * tausqr -(1 / 30) * u * udot * tausqr * tausqr * tau1;
            g1 = tau1 - (1 / 6) * u * tau1 * tausqr - (1 / 12) * udot * tausqr * tausqr -(1 / 120) * u * u * tausqr * tausqr * tau1 -(1 / 120) * u * udot * tausqr * tausqr * tausqr;

            tausqr = tau3 * tau3;

            f3 = 1 - 0.5 * u * tausqr - (1 / 6) * udot * tausqr * tau3 -(1 / 24) * u * u * tausqr * tausqr -(1 / 30) * u * udot * tausqr * tausqr * tau3;
            g3 = tau3 - (1 / 6) * u * tau3 * tausqr - (1 / 12) * udot * tausqr * tausqr -(1 / 120) * u * u * tausqr * tausqr * tau3 -(1 / 120) * u * udot * tausqr * tausqr * tausqr;
        }else {

            theta = angl(&r1, &r2);
            theta1 = angl(&r2, &r3);

            f1 = 1 - ((magr1 * (1 - cos(theta)) / p));
            g1 = (magr1 * magr2 * sin(-theta)) / sqrt(p);
            f3 = 1 - ((magr3 * (1 - cos(theta1)) / p));
            g3 = (magr3 * magr2 * sin(theta1)) / sqrt(p);
        }

        C1 = g3 / (f1 * g3 - f3 * g1);
        C2 = -1;
        C3 = -g1 / (f1 * g3 - f3 * g1);

        double H1 = GM_Earth * tau3 / 12;
        double H3 = -GM_Earth * tau1 / 12;
        double H2 = H1 - H3;

        double G1 = -tau3 / (tau1 * (tau3 - tau1));
        double G3 = -tau1 / (tau3 * (tau3 - tau1));
        double G2 = G1 - G3;

        double D1 = G1 + H1 / pow(magr1, 3);
        double D2 = G2 + H2 / pow(magr2, 3);
        double D3 = G3 + H3 / pow(magr3, 3);

        Matrix C(1,3);
        C(1,1) = C1; C(1,2) = C2; C(1,3) = C3;

        Matrix D(1,3);
        D(1,1) = D1; D(1,2) = D2; D(1,3) = D3;

        double temp = dotProduct(&D, &C)*(-1);
        rhoold1 = temp *(1/(a1 + b1 * u));
        rhoold2 = temp*(-1);
        rhoold3 = temp / (a3 + b3 * u);

        r1 = *Rs1 + Lm1*rhoold1;
        r2 = *Rs2 + Lm2*rhoold2;
        r3 = *Rs3 + Lm3*rhoold3;

    }

    r2 = *Rs1+Lm1*rho1;
    r2 = *Rs2+Lm2*rho2;
    r3 = *Rs3+Lm3*rho3;
}
