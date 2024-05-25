//
// Created by alons on 13/05/2024.
//

#ifndef PROYECTO_GMST_H
#define PROYECTO_GMST_H

#include <cmath>
#include "Frac.h"

//------------------------------------------------------------------------------
// double gmst(double Mjd_UT1)
//------------------------------------------------------------------------------
/**
 * Computes the Greenwich Mean Sidereal Time (GMST) at the given Universal Time
 * (UT1) expressed in Modified Julian Date (MJD).
 *
 * @param Mjd_UT1 Modified Julian Date (UT1).
 *
 * @return Greenwich Mean Sidereal Time (GMST) in radians.
 */
//------------------------------------------------------------------------------
double gmst(double Mjd_UT1);

#endif //PROYECTO_GMST_H
