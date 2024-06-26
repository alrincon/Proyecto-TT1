#include <math.h>
#include <stdio.h>
#include <vector>
#include <iostream>

#include "../include/EFK_Tests.h"

#include "../include/Matrix.h"
#include "../include/R_x.h"
#include "../include/R_y.h"
#include "../include/R_z.h"
#include "../include/AccelPointMass.h"
#include "../include/sign_.h"
#include "../include/AzElPa.h"
#include "../include/Cheb3D.h"
#include "../include/EccAnom.h"
#include "../include/Frac.h"
#include "../include/SAT_const.h"
#include "../include/unit.h"
#include "../include/timediff.h"
#include "../include/Position.h"
#include "../include/Mjday_TDB.h"
#include "../include/Mjday.h"
#include "../include/MeanObliquity.h"
#include "../include/Legendre.h"
#include "../include/IERS.h"
#include "../include/Geodetic.h"
#include "../include/NutAngles.h"
#include "../include/TimeUpdate.h"
#include "../include/NutMatrix.h"
#include "../include/PoleMatrix.h"
#include "../include/PrecMatrix.h"
#include "../include/dotProduct.h"
#include "../include/angl.h"
#include "../include/MeasUpdate.h"
#include "../include/crossProduct.h"
#include "../include/doubler.h"
#include "../include/gibbs.h"
#include "../include/gmst.h"
#include "../include/hgibbs.h"
#include "../include/G_AccelHarmonic.h"
#include "../include/AccelHarmonic.h"

#include "../include/elements.h"
#include "../include/EqnEquinox.h"
#include "../include/G_AccelHarmonic.h"
#include "../include/gast.h"
#include "../include/GHAMatrix.h"
#include "../include/hgibbs.h"
#include "../include/LTC.h"
#include "../include/VarEqn.h"
#include "../include/anglesg.h"
#include "../include/JPL_Eph_DE430.h"
#include "../include/Accel.h"


#define TOL_ 10e-6


using namespace std;


int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) \
  do {                \
    if (!(test)) {    \
      FAIL();         \
      return 1;       \
    }                 \
  } while (0)
#define _verify(test) \
  do {                \
    int r = test();   \
    tests_run++;      \
    if (r) return r;  \
  } while (0)

bool compareVectors (vector<double> v1, vector<double> v2){
    if(v1.size() !=  v2.size()){
        return false;
    }

    for(int i = 0; i < v1.size(); i++){
        if(fabs(v1[i]-v2[i]) >= TOL_){
            return false;
        }
    }
    return true;
}

int R_x_01(){
    double angle = 2.0;
    Matrix sol(3,3);

    sol = R_x(angle);

    _assert(fabs(sol(1,1) -1) < TOL_ && fabs(sol(1,2)) < TOL_ && fabs(sol(1,3)) < TOL_);
    _assert(fabs(sol(2,1)) < TOL_ && fabs(sol(2,2) + 0.416146836547142) < TOL_ && fabs(sol(2,3) - 0.909297426825682) < TOL_);

    return 0;
}

int R_y_01(){
    double angle = 2.0;
    Matrix sol(3,3);

    sol = R_y(angle);

    _assert(fabs(sol(1,1) + 0.416146836547142) < TOL_ && fabs(sol(1,2)) < TOL_ && fabs(sol(1,3) + 0.909297426825682) < TOL_);

    return 0;
}

int R_z_01(){
    double angle = 2.0;
    Matrix sol(3,3);

    sol = R_z(angle);

    _assert(fabs(sol(1,1) + 0.416146836547142) < TOL_ && fabs(sol(1,2) - 0.909297426825682) < TOL_ && fabs(sol(1,3)) < TOL_);

    return 0;
}

int AccelPointMass_01(){
    Matrix r(1, 3);
    r(1,1) = 1;
    r(1,2) = 2;
    r(1,3) = 3;

    Matrix s(1, 3);
    s(1,1) = 0;
    s(1,2) = 4;
    s(1,3) = 6;

    Matrix res = AccelPointMass(r, s, 9.8);

    _assert(compareVectors({res(1,1), res(1,2), res(1,3)}, {-0.18708286933870, 0.26962608631187, 0.40443912946781}));

    return 0;
}

int AzElPa_01() {
    Matrix s(3,1);
    s(1,1) = 1;
    s(2,1) = 2;
    s(3,1) = 3;

    double Az;
    double El;

    Matrix dAds(3,1);
    Matrix dEds(3,1);



    AzElPa(&s, Az, El, dAds, dEds);

    _assert(fabs(Az - 0.463648) < TOL_);
    _assert(fabs(El - 0.930274) < TOL_);
    _assert(compareVectors({dAds(1,1), dAds(2,1), dAds(3,1)}, {0.40000000000000,-0.20000000000000,0.00000000000000}));
    _assert(compareVectors({dEds(1,1), dEds(2,1), dEds(3,1)}, {-0.09583148474999,-0.19166296949998,0.15971914124998}));

    return 0;
}

int sign_01() {
    double a = -2;
    double b = 3;

    double r0 = sign(a,b);
    double r1 = sign(a,-b);

    _assert(fabs(r0 - 2) < TOL_ && fabs(r1 + 2) < TOL_);

    return 0;
}

int Cheb3D_01() {
    double t = 3;
    int N = 5;
    double Ta = 1;
    double Tb = 9;


    Matrix Cx(1,5);
    Cx(1,1) = 1;
    Cx(1,2) = 2;
    Cx(1,3) = 3;
    Cx(1,4) = 4;
    Cx(1,5) = 5;

    Matrix Cy(1,5);
    Cy(1,1) = 1;
    Cy(1,2) = 6;
    Cy(1,3) = 1;
    Cy(1,4) = -4;
    Cy(1,5) = 2;

    Matrix Cz(1,5);
    Cz(1,1) = -1;
    Cz(1,2) = 2;
    Cz(1,3) = 3;
    Cz(1,4) = -1;
    Cz(1,5) = 3;

    Matrix res = Cheb3D(t,N,Ta,Tb,&Cx,&Cy,&Cz);
    _assert(compareVectors({res(1,1), res(1,2), res(1,3)}, {0.00000000000000, -7.50000000000000, -6.00000000000000 }));

    return 0;
}

int EccAnom_01() {
    _assert(abs(EccAnom(2.5,0.7123) - 2.76317) < TOL_);
    _assert(abs(EccAnom(1.76,0.551) - 2.204136) < TOL_);

    return 0;
}

int Frac_01() {
    _assert(abs(Frac(2.56) - 0.56) < TOL_);
    _assert(abs(Frac(0.551) - 0.551) < TOL_);

    return 0;
}

int Geodetic_01() {
    double lon;
    double lat;
    double h;

    Matrix r(3,1);
    r(1,1) = 0.7658;
    r(2,1) = -2.87;
    r(3,1) = 2.2;


    Geodetic(lon, lat, h, &r);

    _assert(compareVectors({lon, lat, h}, {-1.31004, 1.57073, -6356749.41668}));

    return 0;
}

int IERS_01(){
    Matrix eop(13,5);

    eop(1,1) = 0;   eop(1,2) = 0.5;   eop(1,3) = 0;   eop(1,4) = 0;   eop(1,5) = 0;
    eop(2,1) = 3;   eop(2,2) = 1;   eop(2,3) = 7;   eop(2,4) = 0;   eop(2,5) = 0.5;
    eop(3,1) = 4;   eop(3,2) = 0.5;   eop(3,3) = 0;   eop(3,4) = 0.5;   eop(3,5) = 0;
    eop(4,1) = 0.5;   eop(4,2) = 0;   eop(4,3) = 0;   eop(4,4) = 0;   eop(4,5) = 0;
    eop(5,1) = 0;   eop(5,2) = 0.25;   eop(5,3) = -3;   eop(5,4) = 0;   eop(5,5) = 0;
    eop(6,1) = 0;   eop(6,2) = 0;   eop(6,3) = 0;   eop(6,4) = 0;   eop(6,5) = 0;
    eop(7,1) = 0;   eop(7,2) = 1;   eop(7,3) = 10;   eop(7,4) = 0;   eop(7,5) = 0;
    eop(8,1) = 0;   eop(8,2) = -4;   eop(8,3) = 0;   eop(8,4) = 0;   eop(8,5) = 0;
    eop(9,1) = 0;   eop(9,2) = 0;   eop(9,3) = 0;   eop(9,4) = 0;   eop(9,5) = 0;
    eop(10,1) = 0;  eop(10,2) = 23;   eop(10,3) = 0;   eop(10,4) = 0;   eop(10,5) = 0;
    eop(11,1) = 0;   eop(11,2) = -6;   eop(11,3) = 9;   eop(11,4) = 0;   eop(11,5) = 0;
    eop(12,1) = 0;   eop(12,2) = 0;   eop(12,3) = 0.7;   eop(12,5) = 0;   eop(12,5) = 0;
    eop(13,1) = 0;   eop(13,2) = 0.4;   eop(13,3) = 0;   eop(13,5) = 0;   eop(13,5) = 0;

    double Mjd_UTC = 0.6546;
    double x_pole;
    double y_pole;
    double UT1_UTC;
    double LOD;
    double dpsi;
    double deps;
    double dx_pole;
    double dy_pole;
    double TAI_UTC;

    IERS(&eop, Mjd_UTC, 'l', x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
    _assert(compareVectors({x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC}, {-0.00000910213446 ,0.00000000000000 ,6.89140000000000 ,-1.38160000000000 ,0.00000000000000 ,0.00003851456845 ,0.00001851503448 ,0.00000222151325 ,0.40000000000000}));

    IERS(&eop, Mjd_UTC, 'n', x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
    _assert(compareVectors({x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC}, {1.21203e-006,0,1,-4,0,0.000111507,-2.90888e-005,0,0.4}));

    return 0;
}

int Legendre_01(){
    int n = 1;
    int m = 1;
    double fi = 0.3;

    Matrix pnm(n + 1, m + 1);
    Matrix dpnm(n + 1, m + 1);

    Legendre(n, m, fi, pnm, dpnm);

    _assert(compareVectors({pnm(1,1), pnm(1,2), pnm(2,1), pnm(2,2)}, {1.00000000000000, 0.00000000000000, 0.51185601260069, 1.65469133749002}));
    _assert(compareVectors({dpnm(1,1), dpnm(1,2), dpnm(2,1), dpnm(2,2)}, {0.00000000000000, 0.00000000000000, 1.65469133749002, -0.51185601260069}));

    return 0;
}

int MeanObliquity_01(){
    double Mjd_TT = 25.8;

    double res = MeanObliquity(Mjd_TT);
    _assert(abs(res - 0.409413) < TOL_);

    return 0;
}

int Mjday_01(){
    int year = 2017;
    int month = 6;
    int day = 23;
    int hour = 15;
    int min = 45;
    double sec = 25.8;

    double res = Mjday(year, month, day, hour, min, sec);
    _assert(abs(res - 57927.6565486) < TOL_);

    return 0;
}

int Mjday_TBD_01(){
    double Mjd_TT = 0.7689;

    double res  = Mjday_TDB(Mjd_TT);

    _assert(abs(res - 0.7689) < TOL_);
    return 0;
}

int Position_01(){
    double lon = 1.4;
    double lat = 9.9;
    double h = 0.6;


    Matrix res = Position(lon, lat, h);

    _assert(compareVectors({res(1,1), res(2,1), res(3,1)}, {-964624.89547065913212, -5592782.97299871593714, -2900724.09462780551985}));

    return 0;
}

int timediff_01(){
    double UT1_UTC = 0.5;
    double TAI_UTC = 23.7;
    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;

    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

    _assert(compareVectors({UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC}, {-23.2, -4.7, -4.2, 55.884, 4.7}));

    return 0;
}

int unit_01(){
    Matrix v(1,3);
    v(1,1) = 1;
    v(1,2) = -1;
    v(1,3) = 2;

    Matrix res = unit(&v);

    _assert(compareVectors({res(1,1), res(1,2), res(1,3)}, {1/sqrt(6), -1/sqrt(6), 2/sqrt(6)}));

    Matrix v2(1,3);
    v2(1,1) = 0.0000001;
    v2(1,2) = 0.0000001;
    v2(1,3) = 0.0000001;

    Matrix res2 = unit(&v2);

    _assert(compareVectors({res2(1,1), res2(1,2), res2(1,3)}, {0,0,0}));

    return 0;
}

int NutAngles_01(){
    double Mjd_TT = 12.03;
    double dpsi;
    double deps;

    NutAngles(Mjd_TT, dpsi, deps);

    _assert(compareVectors({dpsi, deps}, {3.14455268785733e-05, 3.84461732419543e-05}));

    return 0;
}

int TimeUpdate_01(){
    Matrix P(3,3);
    P(1,1) =  1; P(1,2) = -9; P(1,3) = 10;
    P(2,1) =  4; P(2,2) =  8; P(2,3) =  4;
    P(3,1) =  2; P(3,2) = -5; P(3,3) = -1;


    Matrix Phi(3,3);
    Phi(1,1) =  0; Phi(1,2) =  1; Phi(1,3) =  0;
    Phi(2,1) = -3; Phi(2,2) =  1; Phi(2,3) =  0;
    Phi(3,1) = -1; Phi(3,2) =  0; Phi(3,3) = -1;

    Matrix Qdt(3,3);
    Qdt(1,1) =  1; Qdt(1,2) = -1; Qdt(1,3) =  1;
    Qdt(2,1) =  0; Qdt(2,2) = -1; Qdt(2,3) =  1;
    Qdt(3,1) =  2; Qdt(3,2) =  0; Qdt(3,3) =  1;

    TimeUpdate(P, &Phi, &Qdt);

    _assert(compareVectors({P(1,1), P(1,2), P(1,3), P(2,1), P(2,2), P(2,3), P(3,1), P(3,2), P(3,3)}, {9, -5, -7, 35, 31, 26, 16, 23, 13}));

    TimeUpdate(P, &Phi);
    _assert(compareVectors({P(1,1), P(1,2), P(1,3), P(2,1), P(2,2), P(2,3), P(3,1), P(3,2), P(3,3)}, {31, -74, -61, 46, 22, -55, -18, 57, 31}));

    return 0;
}

int NutMatrix_01(){
    double Mjd_TT = 1.5;
    Matrix res = NutMatrix(Mjd_TT);

    _assert(compareVectors({res(1,1), res(1,2), res(1,3), res(2,1), res(2,2), res(2,3), res(3,1), res(3,2), res(3,3)}, {0.99999999962810, -0.00002501869886, -0.00001085645327, 0.00002501827157, 0.99999999891257, -0.00003935672680, 0.00001085743791, 0.00003935645517, 0.99999999916659}));

    return 0;
}

int PoleMatrix_01(){
    double xp = 1.34;
    double yp = 4.764;

    Matrix res = PoleMatrix(xp, yp);

    _assert(compareVectors({res(1,1), res(1,2), res(1,3), res(2,1), res(2,2), res(2,3), res(3,1), res(3,2), res(3,3)}, {0.22875280780846, -0.97218829537272, 0.05022022759829, 0, 0.05158810997741, 0.99866844693770, -0.97348454169532, -0.22844821130671, 0.01180092500686}));

    return 0;
}

int PrecMatrix_01(){
    double Mjd_1 = 0.2342;
    double Mjd_2 = 0.5452;

    Matrix res = PrecMatrix(Mjd_1, Mjd_2);

    _assert(compareVectors({res(1,1), res(1,2), res(1,3), res(2,1), res(2,2), res(2,3), res(3,1), res(3,2), res(3,3)}, {0.99999999999998, -0.00000019024097, -0.00000008278866, 0.00000019024097, 0.99999999999998, -0.00000000000001, 0.00000008278866, -0.00000000000001, 1}));

    return 0;
}

int dotProduct_01(){
    Matrix v1(1,3);
    v1(1,1) = 1;
    v1(1,2) = 2;
    v1(1,3) = 3;

    Matrix v2(1,3);
    v2(1,1) = 1;
    v2(1,2) = -2;
    v2(1,3) = 0;

    double dot = dotProduct(&v1, &v2);

    _assert(abs(dot + 3) < TOL_);

    return 0;
}

int angl_01(){
    Matrix v1(1,3);
    v1(1,1) = 1;
    v1(1,2) = 2;
    v1(1,3) = 3;

    Matrix v2(1,3);
    v2(1,1) = 1;
    v2(1,2) = -2;
    v2(1,3) = 0;

    double res = angl (&v1, &v2);

    _assert(abs(res - 1.93753) < TOL_);

    return 0;
}

int MeasUpdate_01(){
    Matrix z(3,1);

    z(1,1) =   1;
    z(2,1) =  -1;
    z(3,1) = 0.5;

    Matrix g(3,1);

    g(1,1) = -0.1;
    g(2,1) =  0.7;
    g(3,1) =    3;

    Matrix s(3,1);

    s(1,1) = 1;
    s(2,1) = 0;
    s(3,1) = -1;

    Matrix G(3,3);

    G(1,1) =  1; G(1,2) =  0; G(1,3) =  1;
    G(2,1) =  2; G(2,2) =  1; G(2,3) =  0;
    G(3,1) = -2; G(3,2) = -6; G(3,3) = -1;

    int n = 3;
    Matrix K(3,3);

    K(1,1) =  1; K(1,2) =   3; K(1,3) =  0;
    K(2,1) =  0; K(2,2) =  -7; K(2,3) =  4;
    K(3,1) =  1; K(3,2) =  -2; K(3,3) = -8;

    Matrix x(1,3);

    x(1,1) =  5;
    x(1,2) = -3;
    x(1,3) =  5;

    Matrix P(3,3);

    P(1,1) =  1; P(1,2) = -2; P(1,3) =  5;
    P(2,1) =  1; P(2,2) = -1; P(2,3) =  5;
    P(3,1) =  5; P(3,2) = -1; P(3,3) = -1;


    MeasUpdate(&z, &g, &s, &G , n, K, x, P);

    _assert(compareVectors({K(1,1), K(1,2), K(1,3), K(2,1), K(2,2), K(2,3), K(3,1), K(3,2), K(3,3)},
                            {
        0.08852355358837 ,  0.53588365475814, 0.08883970913690 ,
        -0.17704710717673, -0.07176730951628, -0.17767941827379,
        0.95415744546317 , -0.54536832121404, -0.06386342080304}));

    _assert(compareVectors({x(1,1), x(1,2), x(1,3)},
                           {3.96427442301612, -2.62854884603225, 7.13635788808094}));

    _assert(compareVectors({P(1,1), P(1,2), P(1,3), P(2,1), P(2,2), P(2,3), P(3,1), P(3,2), P(3,3)},
                           {
        0.01612393297502 , -0.03224786595005, 0.07239962061334,
        -0.03224786595005,  0.06449573190009, -0.14479924122668,
        0.08093582042365 , -0.16187164084730, 0.87322162503952
                           }));

    return 0;
}

int crossProduct_01(){
    Matrix v1(1,3);
    v1(1,1) = 1;
    v1(1,2) = 2;
    v1(1,3) = 3;

    Matrix v2(1,3);
    v2(1,1) = 3;
    v2(1,2) = 4;
    v2(1,3) = 5;

    Matrix res = crossProduct(&v1, &v2);

    _assert(compareVectors({res(1,1), res(1,2), res(1,3)}, {-2,4,-2}));


    return 0;
}

int doubler_01(){
    //IN
    double cc1 = 6.5;
    double cc2 = 7.25;
    double magrsite1 = 1.32;
    double magrsite2 = 2.25;
    double magr1in = 0.9;
    double magr2in = 2.8;

    Matrix los1(1, 3);
    los1(1,1) = 1;
    los1(1,2) = -1;
    los1(1,3) = 2;

    Matrix los2(1,3);
    los2(1,1) = 0.2;
    los2(1,2) = 0.1;
    los2(1,3) = 4;

    Matrix los3(1,3);
    los3(1,1) = 1.1;
    los3(1,2) = 2.2;
    los3(1,3) = 3.3;

    Matrix rsite1(1,3);
    rsite1(1,1) = -1;
    rsite1(1,2) = 0;
    rsite1(1,3) = 1;

    Matrix rsite2(1,3);
    rsite2(1,1) = 1;
    rsite2(1,2) = 1;
    rsite2(1,3) = 1;

    Matrix rsite3(1,3);
    rsite3(1,1) = 1;
    rsite3(1,2) = 2.5;
    rsite3(1,3) = 3.6;


    double t1 = 1.5;
    double t3 = 1.8;
    char direct = 'y';

    //OUT
    Matrix r2(1,3);
    Matrix r3(1,3);
    double f1;
    double f2;
    double q1;
    double magr1;
    double magr2;
    double a;
    double deltae32;

    doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, &los1,  &los2,  &los3,  &rsite1,  &rsite2,  &rsite3, t1, t3, direct,  r2, r3,  f1, f2, q1, magr1, magr2, a, deltae32);

    _assert(compareVectors({r2(1,1), r2(1,2), r2(1,3)}, {1.07295049971787, 1.03647524985894, 2.45900999435742}));
    _assert(compareVectors({r3(1,1), r3(1,2), r3(1,3)}, {-0.35732619629596, -0.21465239259192, -0.47197858888789}));

    _assert(abs(f1 - 1.50000003929282) < TOL_);
    _assert(abs(f2 - 1.80000003983532) < TOL_);
    _assert(abs(q1 - 2.34307495852899) < TOL_);
    _assert(abs(magr1 - 1.35488414391266) < TOL_);
    _assert(abs(magr2 - 2.87614913917284) < TOL_);
    _assert(abs(a + 0.33795153370847) < TOL_);
    _assert(abs(deltae32 - 5.05836508736097) < TOL_);

    return 0;
}

int gibbs_01(){
    //IN
    Matrix r1(3,1);
    r1(1,1) = 1;
    r1(2,1) = 2.5;
    r1(3,1) = 3.6;

    Matrix r2(3,1);
    r2(1,1) = 1.5;
    r2(2,1) = 3;
    r2(3,1) = 3.9;

    Matrix r3(3,1);
    r3(1,1) = 2;
    r3(2,1) = 3.3;
    r3(3,1) = 5;

    //OUT
    Matrix v2(3,1);
    double theta;
    double theta1;
    double copa;
    char* error;

    gibbs(&r1, &r2, &r3, v2, theta, theta1, copa, error);

    _assert(compareVectors({v2(1,1), v2(2,1), v2(3,1)}, {450229.09705320186913, 361999.52190559543669, 623055.75882231071591}));

    _assert(abs(theta - 0.08566810351896) < TOL_);
    _assert(abs(theta1 - 0.07374308868996) < TOL_);
    _assert(abs(copa + 0.08373604123886) < TOL_);
    _assert(strcmp (error, "not coplanar") == 0);

    return 0;
}

//CORRECTO
int gmst_01(){
    double Mjd_UT1 = 25.8;
    double gmstime = gmst(Mjd_UT1);

    _assert(abs(gmstime - 0.160403) < TOL_);

    return 0;
}

int hgibbs_01(){
    //IN
    Matrix r1(3,1);
    r1(1,1) = 1;
    r1(2,1) = 2.5;
    r1(3,1) = 3.6;

    Matrix r2(3,1);
    r2(1,1) = 1.5;
    r2(2,1) = 3;
    r2(3,1) = 3.9;

    Matrix r3(3,1);
    r3(1,1) = 2;
    r3(2,1) = 3.3;
    r3(3,1) = 5;

    double Mjd1 = 1.4;
    double Mjd2 = 1.7;
    double Mjd3 = 1.9;

    //OUT
    Matrix v2(3,1);
    double theta;
    double theta1;
    double copa;
    char* error;

    hgibbs(&r1, &r2, &r3, Mjd1, Mjd2, Mjd3, v2, theta, theta1, copa, error);

    _assert(compareVectors({v2(1,1), v2(2,1), v2(3,1)}, {-2645622754722265, -10842213878151156, -13879679200981824}));

    _assert(abs(theta - 0.08566810351896) < TOL_);
    _assert(abs(theta1 - 0.07374308868996) < TOL_);
    _assert(abs(copa + 0.08373604123886) < TOL_);
    _assert(strcmp (error, "angl > 1ø") == 0);

    return 0;
}

int EqnEquinox_01(){
    double Mjd_TT = 12.5;

    double res  = EqnEquinox (Mjd_TT);

    _assert(abs(res - 2.86958e-005) < TOL_);

    return 0;
}

int GHAMatrix_01(){
    double Mjd_UT1 = 12.5;

    Matrix res = GHAMatrix (Mjd_UT1);

    _assert(compareVectors({res(1,1), res(1,2), res(1,3), res(2,1), res(2,2), res(2,3), res(3,1), res(3,2), res(3,3)},
                           {
                                   -0.37326378402851, -0.92772525433596, 0.00000000000000,
                                   0.92772525433596, -0.37326378402851, 0.00000000000000,
                                   0.00000000000000, 0.00000000000000, 1.00000000000000
                           }));

    return 0;
}

int LTC_01(){
    double lon = 11.2;
    double lat = 1.6;

    Matrix res = LTC(lon, lat);

    _assert(compareVectors({res(1,1), res(1,2), res(1,3), res(2,1), res(2,2), res(2,3), res(3,1), res(3,2), res(3,3)},
                           {
                                   0.97917772915132, 0.20300486381875, 0.00000000000000,
                                   -0.20291830316226, 0.97876021074578, -0.02919952230129,
                                   -0.00592764504835, 0.02859152193928, 0.99957360304151
                           }));

    return 0;
}

int elements_01(){
    //IN
    Matrix y(1,6);
    y(1,1) = 1.5;
    y(1,2) = 1.1;
    y(1,3) = 7.2;
    y(1,4) = -1.8;
    y(1,5) = -0.4;
    y(1,6) = 3;

    //OUT
    double p;
    double a;
    double e;
    double i;
    double Omega;
    double omega;
    double M;

    elements (&y, p, a, e, i, Omega, omega, M);

    _assert(abs(p - 8.65399e-013) < TOL_);
    _assert(abs(a - 3.7182) < TOL_);
    _assert(abs(e - 0.999999999999884) < TOL_);
    _assert(abs(i - 1.49643) < TOL_);
    _assert(abs(Omega - 0.340191) < TOL_);
    _assert(abs(omega - 4.47053077244878) < TOL_);
    _assert(abs(M - 3.14159) < TOL_);

    return 0;
}

int gast_01(){
    double Mjd_UT1 = 12.5;

    double res  = gast (Mjd_UT1);

    _assert(abs(res - 4.32986438855874) < TOL_);
    return 0;
}

int VarEqn_01(){
    double x = 0.75;

    Matrix yPhi(1,42);
    yPhi(1,1)= 1;
    yPhi(1,2)= -9;
    yPhi(1,3)= 10;
    yPhi(1,4)= 12;
    yPhi(1,5)= -4;
    yPhi(1,6)= 1;
    yPhi(1,7)= 4;
    yPhi(1,8)= 8;
    yPhi(1,9)= 4;
    yPhi(1,10) = 1;
    yPhi(1,11) = 3;
    yPhi(1,12) = 2;
    yPhi(1,13) = 2;
    yPhi(1,14) = -5;
    yPhi(1,15) = -1;
    yPhi(1,16) = 0;
    yPhi(1,17) = -7;
    yPhi(1,18) = 2;
    yPhi(1,19) = -3;
    yPhi(1,20) = -6;
    yPhi(1,21) = -9;
    yPhi(1,22) = -12;
    yPhi(1,23) = -15;
    yPhi(1,24) = -18;
    yPhi(1,25) = -1;
    yPhi(1,26) = -2;
    yPhi(1,27) = -3;
    yPhi(1,28) = -4;
    yPhi(1,29) = -5;
    yPhi(1,30) = -6;
    yPhi(1,31) = 1.1;
    yPhi(1,32) = 1.2;
    yPhi(1,33) = 1.3;
    yPhi(1,34) = 1.4;
    yPhi(1,35) = 1.5;
    yPhi(1,36) = 1.6;
    yPhi(1,37) = 2.2;
    yPhi(1,38) = 2.3;
    yPhi(1,39) = 2.4;
    yPhi(1,40) = 2.5;
    yPhi(1,41) = 2.6;
    yPhi(1,42) = 2.7;

    Matrix res = VarEqn(0.75, &yPhi);

    _assert(compareVectors({res(1,1), res(1,2), res(1,3), res(1,4), res(1,5), res(1,6), res(1,7), res(1,8), res(1,9), res(1,10), res(1,11), res(1,12), res(1,13), res(1,14), res(1,15), res(1,16), res(1,17), res(1,18), res(1,19), res(1,20), res(1,21), res(1,22), res(1,23), res(1,24), res(1,25), res(1,26), res(1,27), res(1,28), res(1,29), res(1,30), res(1,31), res(1,32), res(1,33), res(1,34), res(1,35), res(1,36), res(1,37), res(1,38), res(1,39), res(1,40), res(1,41), res(1,42)},
                           {12.00000000000000,
                            -4.00000000000000,
                            1.00000000000000,
                            -852697911546758778918763148924580548610343491136233958075266525834128326516284823296584543017407220470099122739514179584.00000000000000,
                            850694796256662385731944628273057283191049735301258602984823586362086918311885479570821767222980146778096147511799447552.00000000000000,
                            -556275986441332890521978914884228922181252922053910468160996699541933594298211398299287070888141738796189351148611698688.00000000000000,
                            1.00000000000000,
                            3.00000000000000,
                            2.00000000000000,
                            -4758975221762135405666726782320935987831503298441363849949608277895693967231641061621385713788689666329719047112464793600.00000000000000,
                            -10331836601268974792181603866044459666545727440049658265511287290062041634639400311820051228970143452544889554368687243264.00000000000000,
                            -10576751299466636314487820714537845140021670039977537079569632555710006201271210813117277079969067836136981980144163880960.00000000000000,
                            0.00000000000000,
                            -7.00000000000000,
                            2.00000000000000,
                            4077203116675019800240667335900375811843652031847065502545268783602715508286507782965322125565342844564504679903983042560.00000000000000,
                            -358283697248133197832519243036920331152183410913399342079756293605711194085301670108626187062555143511005011163981283328.00000000000000,
                            9844422556777667691157343173314111726964377813040726921755302797323930247408914527850381116842613028247310943024622600192.00000000000000,
                            -12.00000000000000,
                            -15.00000000000000,
                            -18.00000000000000,
                            -1884011239681526758841174931272680130374180251345306586029505100259779563920644834472840626579023970566480770065011572736.00000000000000,
                            16406346596350855644993796151017152127636223557862567462534987297097574171403717745129415633297956910827272700015987916800.00000000000000,
                            9836104050456807612293868345879398173157138547338299810955688949557801496633667984975251699784777883986999059624805531648.00000000000000,
                            -4.00000000000000,
                            -5.00000000000000,
                            -6.00000000000000,
                            -628003746560508728488840311251600742556318808607857128962752602729832876507811662580707201485940651020614164177863311360.00000000000000,
                            5468782198783618548331265383672384042545407852620855820844995765699191390467905915043138544432652303609090900005329305600.00000000000000,
                            3278701350152269395222840781132425358620787457620011670032312080542694143676959277235323573968659967163879112385409056768.00000000000000,
                            1.40000000000000,
                            1.50000000000000,
                            1.60000000000000,
                            -106575424611410456472108652614108526587074418481326670401189109893304459138300324000622023185840913189989888784152395776.00000000000000,
                            -3076626255674850926640509830935634780679349933336955669639598866794913138448450690925565404955052103130485049688825266176.00000000000000,
                            -1396388191150489424652511708825047380680523475128393103242748453978477767342075592185105076161978088793633116550571491328.00000000000000,
                            2.50000000000000,
                            2.60000000000000,
                            2.70000000000000,
                            -292888803805617965311891242025631567108850783581291276958626259556931306526063263332054961401067053233498184555713003520.00000000000000,
                            -5859349095050989708353147049285899858533953351949742044747954294940799845946370585711948741487700487683512193889937129472.00000000000000,
                            -2571758052899277745466011224395160620033051523055774889724519374204504682847902133446717370876335937475263282571181555712.00000000000000}));

    return 0;
}

int AccelHarmonic_01(){
    int n_max = 5;
    int m_max = 6;

    Matrix E(3,3);

    E(1,1) =  1; E(1,2) =   3; E(1,3) =  0;
    E(2,1) =  0; E(2,2) =  -7; E(2,3) =  4;
    E(3,1) =  1; E(3,2) =  -2; E(3,3) = -8;

    Matrix r(3,1);

    r(1,1) = 1;
    r(2,1) = 1.5;
    r(3,1) = 2.3;

    Matrix res = AccelHarmonic(&r, &E, n_max, m_max);

    _assert(compareVectors({res(1,1), res(2,1), res(3,1)}, {-22463029536785409205627661515227136.00000000000000,34680465296512803264427921359503360.00000000000000,111788941618518297131034606995243008.00000000000000}));


    return 0;
}

int G_AccelHarmonic_01(){
    int n_max = 5;
    int m_max = 6;

    Matrix r(3,1);
    r(1,1) = 1;
    r(2,1) = 1.5;
    r(3,1) = 2.3;

    Matrix U(3,3);

    U(1,1) =  1; U(1,2) =   3; U(1,3) =  0;
    U(2,1) =  0; U(2,2) =  -7; U(2,3) =  4;
    U(3,1) =  1; U(3,2) =  -2; U(3,3) = -8;

    Matrix res = G_AccelHarmonic(&r, &U, n_max, m_max);

    _assert(compareVectors({res(1,1), res(1,2), res(1,3), res(2,1), res(2,2), res(2,3), res(3,1), res(3,2), res(3,3)},
                           {
                                   -7505249416737317816258887975698432.00000000000000, 9552330075736158590841992498905088.00000000000000, 65669745033203205579116791268900864.00000000000000,
                                   21827074621157776549830003873808384.00000000000000, 164449896652331264655708328759918592.00000000000000, -340715861701754549157210223089287168.00000000000000,
                                   57464944577628769597750538693246976.00000000000000, -124948575499782943223932591582216192.00000000000000, -164786823477902409639772383447351296.00000000000000
                           }));

    return 0;
}

int anglessg_01(){
    double az1 = 1.05591;
    double az2 = 1.3631;
    double az3 = 1.97616;
    double el1 = 0.282625;
    double el2 = 0.453435;
    double el3 = 0.586427;
    double Mjd1 = 49746.1;
    double Mjd2 = 49773.1;
    double Mjd3 = 49782.1;

    Matrix Rs(3,1);
    Rs(1,1) = -5512567.84003606811166;
    Rs(2,1) = -2196994.44666933314875;
    Rs(3,1) = 2330804.96614688728005;

    Matrix r2(3,1);
    Matrix v2(3,1);

    anglesg(az1, az2, az3, el1, el2, el3,Mjd1,Mjd2,Mjd3,&Rs,&Rs,&Rs, r2, v2);

    _assert(compareVectors({r2(1,1), r2(2,1), r2(3,1)}, {7804856167655.16113281250000,-39421404795886.11718750000000,-14211923962295.30664062500000}));
    _assert(compareVectors({v2(1,1),v2(2,1), v2(3,1)}, {3.63245858917435,-18.34534848119620,-6.61300506702719}));


    return 0;
}

int JPL_Eph_DE430(){
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

    JPL_Eph_DE430(49746.1128763999,r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, r_Neptune,r_Pluto,r_Moon,r_Sun);

    _assert(compareVectors({r_Mercury(1,1),r_Mercury(1,2), r_Mercury(1,3)}, {83772905189.86315917968750, -65290454439.25350189208984, -23390629645.41803741455078}));

    _assert(compareVectors({r_Venus(1,1),r_Venus(1,2), r_Venus(1,3)}, {-15227808349.19389343261719,-110136124950.86943054199219,-41022337370.63686370849609}));

    _assert(compareVectors({r_Earth(1,1),r_Earth(1,2), r_Earth(1,3)}, {-92472772680.09974670410156,106393603655.53057861328125,46129569866.12557983398438}));

    _assert(compareVectors({r_Mars(1,1),r_Mars(1,2), r_Mars(1,3)}, {-88277798473.76695251464844,46964990228.49914550781250,29071127006.03056335449219}));

    _assert(compareVectors({r_Jupiter(1,1),r_Jupiter(1,2), r_Jupiter(1,3)}, {-298383266667.87084960937500,-754497339800.52270507812500,-314410139017.64038085937500}));

    _assert(compareVectors({r_Saturn(1,1),r_Saturn(1,2), r_Saturn(1,3)}, {1482035974927.73242187500000,-453870923240.66571044921875,-249401413783.27105712890625}));

    _assert(compareVectors({r_Uranus(1,1),r_Uranus(1,2), r_Uranus(1,3)}, {1412370256343.40747070312500,-2511353537383.14208984375000,-1118108003556.99047851562500}));

    _assert(compareVectors({r_Neptune(1,1),r_Neptune(1,2), r_Neptune(1,3)}, {1871252960132.36352539062500,-3928974841842.59619140625000,-1655019851828.05908203125000}));

    _assert(compareVectors({r_Pluto(1,1),r_Pluto(1,2), r_Pluto(1,3)}, {-2171412614553.17797851562500,-3915432033426.32910156250000,-552715859865.25683593750000}));

    _assert(compareVectors({r_Moon(1,1),r_Moon(1,2), r_Moon(1,3)}, {89463061.41452670097351,-336587597.53726130723953,-114638183.88342174887657}));

    _assert(compareVectors({r_Sun(1,1),r_Sun(1,2), r_Sun(1,3)}, {92300062279.84397888183594,-105373881161.06535339355469,-45685797078.33899688720703}));


    return 0;
}

int Accel_01(){
    double x = 0.5;
    Matrix Y(1,6);
    Y(1,1) = 1;
    Y(1,2) = 2;
    Y(1,3) = 3;
    Y(1,4) = 4;
    Y(1,5) = 5;
    Y(1,6) = 6;

    Matrix res = Accel(x, &Y);

    _assert(compareVectors({res(1,1),res(1,2), res(1,3), res(1,4), res(1,5), res(1,6)}, {
    4.00000000000000,
     5.00000000000000,
     6.00000000000000,
     -471010558047621988587913909448671909650375564787349795468405082677115760841073283365019401699788064431255484852499948353752305500160.00000000000000,
    -2246660984646394166059964940525112324986811344270867184440307547912545646199989268019447929414451705624083443865098096881241719570432.00000000000000,
    -3731402515517448595593712178093056039286968131588511044842860459458435608674245892077054060650874330524857460759519957201011613892608.00000000000000}));


    return 0;
}

int all_tests(){
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);
    _verify(AccelPointMass_01);
    _verify(sign_01);
    _verify(AzElPa_01);
    _verify(Cheb3D_01);
    _verify(EccAnom_01);
    _verify(Frac_01);
    _verify(unit_01);
    _verify(timediff_01);
    _verify(Position_01);
    _verify(Mjday_TBD_01);
    _verify(Mjday_01);
    _verify(MeanObliquity_01);
    _verify(Legendre_01);
    _verify(IERS_01);
    _verify(Geodetic_01);
    _verify(NutAngles_01);
    _verify(TimeUpdate_01);
    _verify(NutMatrix_01);
    _verify(PoleMatrix_01);
    _verify(PrecMatrix_01);
    _verify(dotProduct_01);
    _verify(angl_01);
    _verify(MeasUpdate_01);
    _verify(crossProduct_01);
    _verify(doubler_01);
    _verify(gibbs_01);
    _verify(gmst_01);
    _verify(hgibbs_01);
    _verify(EqnEquinox_01);
    _verify(GHAMatrix_01);
    _verify(LTC_01);
    _verify(elements_01);
    _verify(gast_01);
    _verify(AccelHarmonic_01);
    _verify(G_AccelHarmonic_01);
    _verify(VarEqn_01);
    _verify(anglessg_01);
    _verify(JPL_Eph_DE430);
    _verify(Accel_01);

    return 0;
}

extern Matrix eopdata;
extern Matrix obs;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix PC;
extern aux AuxParam;


int main(){
    int nobs = 46;
    int infFile = 21413;

    //fichero DE430Coeff.txt
    FILE *f1 = fopen("../data/DE430Coeff.txt", "r");
    if (!f1) {
        std::cerr << "No se pudo abrir el archivo." << std::endl;
        return 1;
    }

    for (int i = 1; i <= 2285; ++i) {
        for (int j = 1; j <= 1020; ++j) {
            if (fscanf(f1, "%lf", &PC(i, j)) != 1) {
                std::cerr << "Error al leer el archivo." << std::endl;
                fclose(f1);
                return 1;
            }
        }
    }

    fclose(f1);

    //fichero GGM03S.txt
    FILE *fp2 = fopen("../data/GGM03S.txt", "r");
    if (fp2 == nullptr) {
        printf("error");
        exit(1);
    }

    for (int i = 0; i < 181; i++) {
        int valor1;
        int valor2;
        double valor5;
        double valor6;

        for(int j = 0; j <= i; j++){
            fscanf(fp2, "%d %d %lf %lf %lf %lf", &valor1, &valor2, &Cnm(i+1, j+1), &Snm(i+1, j+1), &valor5, &valor6);
        }
    }

    if (fclose(fp2))   // Close the stream.
        perror("fclose error");
    else printf("File mylib/myfile closed successfully.\n");

    //parametros de orientación de la tierra
    FILE *fp3 = fopen("../data/eop19620101.txt", "r");
    if (fp3 == nullptr) {
        printf("error");
        exit(1);
    }

    for (int i = 1; i <= infFile; i++) {
        fscanf(fp3, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &eopdata( 1,i), &eopdata( 2,i), &eopdata( 3,i),
               &eopdata(4,i), &eopdata(5,i), &eopdata(6,i),
               &eopdata(7,i), &eopdata(8,i), &eopdata(9,i),
               &eopdata(10,i), &eopdata(11,i), &eopdata(12,i),
               &eopdata(13,i));
    }
    if (fclose(fp3))   // Close the stream.
        perror("fclose error");
    else printf("File mylib/myfile closed successfully.\n");


    //observaciones
    FILE *fp4 = fopen("../data/GEOS3.txt", "r");
    if (fp4 == nullptr) {
        printf("error");
        exit(1);
    }
    int year;
    int month;
    int day;

    int hour;
    int min;
    double sec;

    double az;
    double el;
    double Dist;

    int scancount = 1;
    while (fscanf(fp4, "%d/%d/%d %d:%d:%lf %lf %lf %lf",&year, &month, &day, &hour, &min, &sec,&az, &el, &Dist) == 9 && scancount <= nobs) {
        obs(scancount,1) = Mjday(year, month, day, hour,min, sec);
        obs(scancount,2) = Rad*az;
        obs(scancount,3) = Rad*el;
        obs(scancount,4) = 1e3*Dist;
        scancount++;
    }

    if (fclose(fp4))   // Close the stream.
        perror("fclose error");
    else printf("File mylib/myfile closed successfully.\n");

    double Mjd_UTC = obs(9,1);

    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_UTC;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;
    //parametros de modelo

    int result = all_tests();
    if(result == 0){
        printf("PASSED\n");
    }

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}