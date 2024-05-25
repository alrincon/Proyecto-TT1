#include "../include/realmod.h"

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
double realmod(double x, double y)
{
    if (y == 0 )
        return x;
    else if ( x == y)
        return 0.0;
    else
        return x - floor(x/y) *y;
}