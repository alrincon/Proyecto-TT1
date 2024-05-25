#ifndef PROYECTO_REALMOD_H
#define PROYECTO_REALMOD_H
#include <cmath>

//------------------------------------------------------------------------------
// double realmod(double x, double y)
//------------------------------------------------------------------------------
/**
 * Computes the remainder of x divided by y, similar to the modulus but ensures
 * the result is always non-negative.
 *
 * @param x The dividend.
 * @param y The divisor.
 * @return The remainder of x divided by y. If y is zero, returns x.
 */
//------------------------------------------------------------------------------
double realmod(double x, double y);

#endif
