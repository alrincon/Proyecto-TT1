#include "../include/Frac.h"

//------------------------------------------------------------------------------
// Frac(double x)
//------------------------------------------------------------------------------
/**
 * Computes the fractional part of a number.
 *
 * This function computes the fractional part of a given number.
 *
 * @param x Input number.
 * @return Fractional part of the input number.
 */
//------------------------------------------------------------------------------
double Frac(double x) {
    return x - floor(x);
}