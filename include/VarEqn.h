#ifndef PROYECTO_VAREQN_H
#define PROYECTO_VAREQN_H

#include "Matrix.h"
#include "SAT_const.h"
#include "IERS.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "PoleMatrix.h"
#include "NutMatrix.h"
#include "GHAMatrix.h"
#include "AccelHarmonic.h"
#include "G_AccelHarmonic.h"
#include "types.h"

Matrix VarEqn(double x, Matrix* yPhi);

#endif
