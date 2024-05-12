#ifndef PROYECTO_GIBBS_H
#define PROYECTO_GIBBS_H

#include "Matrix.h"
#include "SAT_const.h"
#include "unit.h"
#include "crossProduct.h"
#include "dotProduct.h"
#include "angl.h"
#include <stdio.h>
#include <string.h>

void gibbs(Matrix *r1, Matrix *r2, Matrix *r3, Matrix &v2, double &theta, double &theta1, double &copa, char* &error);

#endif
