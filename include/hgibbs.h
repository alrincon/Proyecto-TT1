#ifndef PROYECTO_HGIBBS_H
#define PROYECTO_HGIBBS_H

#include "Matrix.h"
#include "SAT_const.h"
#include <stdio.h>
#include <string.h>
#include "crossProduct.h"
#include "dotProduct.h"
#include "unit.h"

void hgibbs (Matrix* r1, Matrix* r2, Matrix* r3, double Mjd1, double Mjd2, double Mjd3, Matrix& v2, double &theta, double &theta1,double &copa, char* &error);

#endif
