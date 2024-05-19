#include "../include/anglesdr.h"


void anglesdr (double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, Matrix *rsite1, Matrix *rsite2, Matrix *rsite3, Matrix &r2, Matrix &v2){
    double magr1in = 1.1*R_Earth;
    double magr2in = 1.11*R_Earth;

    char* direct = static_cast<char *>(malloc(100));

    strcpy(direct,"y");

    double tol    = 1e-8*R_Earth;
    double pctchg = 0.005;

    double t1 = (Mjd1 - Mjd2)*86400.0;
    double t3 = (Mjd3 - Mjd2)*86400.0;

    Matrix los1(1,3);
    los1(1,1) = cos(el1)*sin(az1);
    los1(1,2) = cos(el1)*cos(az1);
    los1(1,3) = sin(el1);

    Matrix los2(1,3);
    los2(1,1) = cos(el2)*sin(az2);
    los2(1,2) = cos(el2)*cos(az2);
    los2(1,3) = sin(el2);

    Matrix los3(1,3);
    los3(1,1) = cos(el3)*sin(az3);
    los3(1,2) = cos(el3)*cos(az3);
    los3(1,3) = sin(el3);

    double lon1, lat1, h1;
    Geodetic(lon1, lat1, h1, rsite1);

    double lon2, lat2, h2;
    Geodetic(lon2, lat2, h2, rsite2);

    double lon3, lat3, h3;
    Geodetic(lon3, lat3, h3, rsite3);

    Matrix M1 = LTC(lon1, lat1);
    Matrix M2 = LTC(lon2, lat2);
    Matrix M3 = LTC(lon3, lat3);

    // body-fixed system
    los1 = M1.transpose()*los1;
    los2 = M1.transpose()*los2;
    los3 = M1.transpose()*los3;

    //mean of date system (J2000)
    double Mjd_UTC = Mjd1;
    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;
    IERS(Global::eopdata,Mjd_UTC,'l', x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);

    double UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;
    timediff(UT1_UTC,TAI_UTC, UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);

    double Mjd_TT = Mjd_UTC + TT_UTC/86400;
    double Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    Matrix P = PrecMatrix(MJD_J2000,Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    los1 = E.transpose()*los1;
    Matrix t = E.transpose()*(*rsite1);
    rsite1 = &t;

    //2

    Mjd_UTC = Mjd2;
    IERS(Global::eopdata,Mjd_UTC,'l', x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);

    timediff(UT1_UTC,TAI_UTC, UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);

    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    P = PrecMatrix(MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    los2 = E.transpose()*los2;
    t = E.transpose()*(*rsite2);
    rsite2 = &t;

    //3

    Mjd_UTC = Mjd3;
    IERS(Global::eopdata,Mjd_UTC,'l', x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);

    timediff(UT1_UTC,TAI_UTC, UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);

    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    P = PrecMatrix(MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    los3 = E.transpose()*los3;
    t = E.transpose()*(*rsite3);
    rsite3 = &t;

    //
    double magr1old  = 99999999.9;
    double magr2old  = 99999999.9;
    double magrsite1 = rsite1->norm();
    double magrsite2 = rsite2->norm();
    double magrsite3 = rsite3->norm();

    double cc1 = dotProduct(&los1,rsite1)*2.0;
    double cc2 = 2.0*dotProduct(&los2,rsite2);
    double ktr = 0;

    while (abs(magr1in-magr1old) > tol || abs(magr2in-magr2old) > tol) {
        ktr = ktr + 1;
        Matrix r2(1,3);
        Matrix r3(1,3);
        double f1, f2, q1, magr1, magr2, a, deltae32;
        doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, &los1, &los2, &los3, rsite1, rsite2, rsite3, t1, t3, *direct, r2, r3, f1, f2, q1, magr1, magr2, a, deltae32);

        double f = 1.0 - a / magr2 * (1.0 - cos(deltae32));
        double g = t3 - sqrt(pow(a, 3) /GM_Earth)*(deltae32 - sin(deltae32));
        Matrix v2 = (r3 - r2*f) * (1/g);

        double magr1o = magr1in;
        magr1in = (1.0 + pctchg) * magr1in;
        double deltar1 = pctchg * magr1in;

        double f1delr1,f2delr1,q2;
        doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, &los1, &los2, &los3, rsite1, rsite2, rsite3, t1, t3, *direct, r2,r3,f1delr1,f2delr1,q2,magr1,magr2,a,deltae32);
        double pf1pr1 = (f1delr1 - f1) / deltar1;
        double pf2pr1 = (f2delr1 - f2) / deltar1;

        magr1in = magr1o;
        deltar1 = pctchg * magr1in;
        double magr2o = magr2in;
        magr2in = (1.0 + pctchg) * magr2in;
        double deltar2 = pctchg * magr2in;

        double f1delr2, f2delr2, q3;

        doubler(cc1, cc2, magrsite1, magrsite2, magr1in,magr2in, &los1, &los2, &los3, rsite1, rsite2, rsite3, t1, t3, *direct, r2, r3, f1delr2, f2delr2, q3, magr1, magr2, a, deltae32);
        double pf1pr2 = (f1delr2 - f1) / deltar2;
        double pf2pr2 = (f2delr2 - f2) / deltar2;

        magr2in = magr2o;
        deltar2 = pctchg * magr2in;

        double delta = pf1pr1 * pf2pr2 - pf2pr1 * pf1pr2;
        double delta1 = pf2pr2 * f1 - pf1pr2 * f2;
        double delta2 = pf1pr1 * f2 - pf2pr1 * f1;

        deltar1 = -delta1 / delta;
        deltar2 = -delta2 / delta;

        magr1old = magr1in;
        magr2old = magr2in;

        magr1in = magr1in + deltar1;
        magr2in = magr2in + deltar2;

    }

    Matrix r3(1,3);
    double f1,f2,q1,magr1,magr2,a,deltae32;
    doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,&los1,&los2,&los3,rsite1,rsite2,rsite3,t1,t3,*direct,r2,r3,f1,f2,q1,magr1,magr2,a,deltae32);

    double f  = 1.0 - a/magr2*(1.0-cos(deltae32));
    double g  = t3 - sqrt(pow(a,3)/GM_Earth)*(deltae32-sin(deltae32));
    v2 = (r3 - r2*f)*(1/g);
}