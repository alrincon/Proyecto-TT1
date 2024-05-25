#ifndef PROYECTO_NUTANGLES_H
#define PROYECTO_NUTANGLES_H


#include <iostream>
#include <vector>
#include "SAT_const.h"
#include <cmath>

//------------------------------------------------------------------------------
// void NutAngles(double Mjd_TT, double& dpsi, double& deps)
//------------------------------------------------------------------------------
/**
 * Computes the nutation in longitude and obliquity.
 *
 * @param Mjd_TT  Modified Julian Date (MJD) in Terrestrial Time (TT).
 * @param dpsi    Nutation in longitude [rad].
 * @param deps    Nutation in obliquity [rad].
 */
//------------------------------------------------------------------------------
void  NutAngles (double Mjd_TT, double& dpsi, double& deps);


#endif //PROYECTO_NUTANGLES_H
