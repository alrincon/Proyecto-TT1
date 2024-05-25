#ifndef PROYECTO_GAST_H
#define PROYECTO_GAST_H

#include <cmath>
#include "EqnEquinox.h"
#include "gmst.h"
#include "realmod.h"

//------------------------------------------------------------------------------
// gast(double Mjd_UT1)
//------------------------------------------------------------------------------
/**
 * Computes the Greenwich Apparent Sidereal Time (GAST).
 *
 * This function calculates the Greenwich Apparent Sidereal Time (GAST) given
 * the Modified Julian Date in the UT1 time scale.
 *
 * @param Mjd_UT1 Modified Julian Date in the UT1 time scale.
 * @return Greenwich Apparent Sidereal Time (GAST).
 */
//------------------------------------------------------------------------------
double gast (double Mjd_UT1);

#endif
