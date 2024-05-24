#ifndef PROYECTO_DEINTEG_H
#define PROYECTO_DEINTEG_H

#include "Matrix.h"
#include "sign_.h"

Matrix DEInteg(Matrix (*func)(double, Matrix*), double t, double tout, double relerr, double abserr, int n_eqn, Matrix *y);
#endif
