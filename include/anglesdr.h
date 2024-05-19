#ifndef PROYECTO_ANGLESDR_H
#define PROYECTO_ANGLESDR_H

#include "Matrix.h"
#include "SAT_const.h"
#include "Global.h"
#include <stdio.h>
#include <string.h>
#include "Geodetic.h"
#include "LTC.h"
#include "IERS.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "dotProduct.h"
#include "doubler.h"


void anglesdr (double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, Matrix *rsite1, Matrix *rsite2, Matrix *rsite3, Matrix &r2, Matrix &v2);

#endif
