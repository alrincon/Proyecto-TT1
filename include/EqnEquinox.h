#ifndef PROYECTO_EQNEQUINOX_H
#define PROYECTO_EQNEQUINOX_H

#include <cmath>
#include "NutAngles.h"
#include "MeanObliquity.h"

//------------------------------------------------------------------------------
// EqnEquinox(double Mjd_TT)
//------------------------------------------------------------------------------
/**
 * Computes the equation of the equinoxes.
 *
 * This function computes the equation of the equinoxes, which is the
 * difference between apparent and mean sidereal times, using the nutation in
 * longitude and the mean obliquity of the ecliptic.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time).
 * @return Equation of the equinoxes.
 */
//------------------------------------------------------------------------------
double EqnEquinox (double Mjd_TT);

#endif
