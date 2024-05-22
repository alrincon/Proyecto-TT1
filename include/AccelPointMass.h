#ifndef ACCELPOINTMASS_H
#define ACCELPOINTMASS_H

#include "../include/Matrix.h"
#include <cmath>
using namespace std;

Matrix AccelPointMass(Matrix& r, Matrix& s, double GM);
Matrix AccelPointMassT(Matrix& r, Matrix& s, double GM);

#endif
