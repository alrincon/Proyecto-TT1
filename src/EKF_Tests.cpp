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

    return 0;
}



int main(){
    int result = all_tests();
    cout<<SAT_Const::R_Earth<<endl;
    if(result == 0){
        printf("PASSED\n");
    }

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
