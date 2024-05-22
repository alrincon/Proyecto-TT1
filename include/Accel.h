#ifndef PROYECTO_ACCEL_H
#define PROYECTO_ACCEL_H

#include "Matrix.h"
#include "IERS.h"
#include "timediff.h"
#include "SAT_const.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "Mjday_TDB.h"
#include "AccelHarmonic.h"
#include "AccelPointMass.h"
#include "JPL_Eph_DE430.h"
#include "types.h"


Matrix Accel(double x, Matrix* Y);

Matrix AccelT(double x, Matrix* Y);
#endif
