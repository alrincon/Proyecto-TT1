#include <math.h>
#include <stdio.h>
#include <vector>
#include <iostream>

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
    double s[3] = {1, 2, 3};
    double Az, El;
    double dAds[3], dEds[3];

    AzElPa(s, Az, El, dAds, dEds);

    _assert(fabs(Az - 0.463647609000806) < TOL_ && fabs(El-0.930274014115472)<TOL_);

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
    Matrix Cx(1,5,new double[5]{0.4,0.1,0.5,0.6,0.7},5);
    Matrix Cy(1,5,new double[5]{-0.2,-0.3,-0.5,-0.1,0.1},5);
    Matrix Cz(1,5, new double[5]{0.9,1.1,1.5,1.1,0.7},5);

    Matrix res = Cheb3D(t,N,Ta,Tb,Cx,Cy,Cz);
    _assert(compareVectors({res(1,1), res(1,2), res(1,3)}, {0.3500, 0.0500, 0.3500}));

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

//checkear que los valores son correctos
int Geodetic_01() {
    double lon;
    double lat;
    double h;

    Matrix r(1,3);
    r(1,1) = 0.7658;
    r(1,2) = -2.87;
    r(1,3) = 2.2;

    Geodetic(lon, lat, h, &r);

    _assert(compareVectors({lon, lat, h}, {-1.31004, 1.57073, -6356749.41668}));

    return 0;
}

//checkear que los valores son correctos
int IERS_01(){
    Matrix eop(13,5);

    eop(1,1) = 0;   eop(1,2) = 0.5;   eop(1,3) = 0;   eop(1,4) = 0;   eop(1,5) = 0;
    eop(2,1) = 3;   eop(2,2) = 1;   eop(2,3) = 7;   eop(2,4) = 0;   eop(2,5) = 0.5;
    eop(3,1) = 4;   eop(3,2) = 0.5;   eop(3,3) = 0;   eop(3,4) = 0.5;   eop(3,5) = 0;
    eop(4,1) = 0.5;   eop(4,2) = 0;   eop(4,3) = 0;   eop(4,4) = 0;   eop(4,5) = 0;
    eop(5,1) = 0;   eop(5,2) = 0;   eop(5,3) = 0;   eop(5,4) = 0;   eop(5,5) = 0;
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

    IERS(&eop, Mjd_UTC, 'n', x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
    _assert(compareVectors({x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC}, {0.00000121203420, 0, 1, -4, 0, 0.00011150714666, -0.00002908882087, 0, 0.4}));

    IERS(&eop, Mjd_UTC, 'l', x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
    _assert(compareVectors({x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC}, {-0.00000910213446, 0, 6.89140000000000, -1.38160000000000, 0, 0.00003851456845, 0.00001851503448, 0.00000222151325, 0.4}));

    return 0;
}

//checkear si los valores son correctos
int Legendre_01(){
    int n = 1;
    int m = 1;
    double fi = 0.3;

    Matrix pnm(n + 1, m + 1);
    Matrix dpnm(n + 1, m + 1);

    Legendre(n, m, fi, pnm, dpnm);

    _assert(compareVectors({pnm(1,1), pnm(1,2), pnm(2,1), pnm(2,2)}, {1, 0, 0, 1.65469133749002}));
    _assert(compareVectors({dpnm(1,1), dpnm(1,2), dpnm(2,1), dpnm(2,2)}, {0, 0, 0, -0.51185601260069}));

    return 0;
}

//checkear si los valores son correctos
int MeanObliquity_01(){
    double Mjd_TT = 25.8;

    double res = MeanObliquity(Mjd_TT);
    _assert(abs(res - 0.409413) < TOL_);

    return 0;
}

//checkear si los valores son correctos
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

//checkear si los valores son correctos
int Mjday_TBD_01(){
    double Mjd_TT = 0.7689;

    double res  = Mjday_TDB(Mjd_TT);

    _assert(abs(res - 0.7689) < TOL_);
    return 0;
}

//checkear si los valores son correctos
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

//checkear si los valores son correctos
int NutAngles_01(){
    double Mjd_TT = 12.03;
    double dpsi;
    double deps;

    NutAngles(Mjd_TT, dpsi, deps);
    _assert(compareVectors({dpsi, deps}, {-6.35126e-011, 7.6151e-011}));

    return 0;
}

//checkear si los valores son correctos
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

//checkear si los valores son correctos
int NutMatrix_01(){
    double Mjd_TT = 1.5;
    Matrix res = NutMatrix(Mjd_TT);

    _assert(compareVectors({res(1,1), res(1,2), res(1,3), res(2,1), res(2,2), res(2,3), res(3,1), res(3,2), res(3,3)}, {1, -0.00000000102292, -0.00000000044388, 0.00000000102292, 1, -0.00000000077570, 0.00000000044388, 0.00000000077570, 1}));

    return 0;
}

//checkear si los valores son correctos
int PoleMatrix_01(){
    double xp = 1.34;
    double yp = 4.764;

    Matrix res = PoleMatrix(xp, yp);

    _assert(compareVectors({res(1,1), res(1,2), res(1,3), res(2,1), res(2,2), res(2,3), res(3,1), res(3,2), res(3,3)}, {0.22875280780846, -0.97218829537272, 0.05022022759829, 0, 0.05158810997741, 0.99866844693770, -0.97348454169532, -0.22844821130671, 0.01180092500686}));

    return 0;
}

//checkear si los valores son correctos
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

//checkear si los valores son correctos
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

//checkear si los valores son correctos
int MeasUpdate_01(){
    Matrix z(1,3);

    z(1,1) =   1;
    z(1,2) =  -1;
    z(1,3) = 0.5;

    Matrix g(1,3);

    g(1,1) = -0.1;
    g(1,2) =  0.7;
    g(1,3) =    3;

    Matrix s(1,3);

    s(1,1) = 1;
    s(1,2) = 0;
    s(1,3) = -1;

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
        0.95415744546317 , -0.54536832121404, -0.06386342080304
                            }));

    _assert(compareVectors({x(1,1), x(1,2), x(1,3)},
                           {5.08852355358837, -2.46411634524186, 5.08883970913690}));

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

//checkear si los valores son correctos
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
    _assert(abs(a + 0.05473104705247) < TOL_);
    _assert(abs(deltae32 + 6.78261898699634) < TOL_);

    return 0;
}

//checkear si los valores son correctos
int gibbs_01(){
    //IN
    Matrix r1(1,3);
    r1(1,1) = 1;
    r1(1,2) = 2.5;
    r1(1,3) = 3.6;

    Matrix r2(1,3);
    r2(1,1) = 1.5;
    r2(1,2) = 3;
    r2(1,3) = 3.9;

    Matrix r3(1,3);
    r3(1,1) = 2;
    r3(1,2) = 3.3;
    r3(1,3) = 5;

    //OUT
    Matrix v2(1,3);
    double theta;
    double theta1;
    double copa;
    char* error;

    gibbs(&r1, &r2, &r3, v2, theta, theta1, copa, error);

    _assert(compareVectors({v2(1,1), v2(1,2), v2(1,3)}, {450229.09705320186913, 361999.52190559543669, 623055.75882231071591}));

    _assert(abs(theta - 0.08566810351896) < TOL_);
    _assert(abs(theta1 - 0.07374308868996) < TOL_);
    _assert(abs(copa + 0.08373604123886) < TOL_);
    _assert(strcmp (error, "not coplanar") == 0);

    return 0;
}

//checkear si los valores son correctos
int gmst_01(){
    double Mjd_UT1 = 25.8;
    double gmstime = gmst(Mjd_UT1);

    _assert(abs(gmstime - 0.160403) < TOL_);

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

    return 0;
}



int main(){
    int result = all_tests();
    if(result == 0){
        printf("PASSED\n");
    }

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
