#include <math.h>
#include <stdio.h>

#include "../include/Matrix.h"
#include "../include/R_z.h"
#include "string.h"

#define TOL_ 10e-14


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


int R_z_01(){
    double angle = 2.0;
    Matrix sol(3,3);

    sol = R_z(angle);
    sol = R_z(angle);
    sol.print();

    _assert(fabs(sol(1,1)) -0.416146836547142 < TOL_ && fabs(sol(1,2)) - 0.909297426825682 < TOL_ && fabs(sol(1,3)) < TOL_);
}

int main(){
    _verify(R_z_01);
}