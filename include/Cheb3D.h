#ifndef CHEB3D_H
#define CHEB3D_H

#include "Matrix.h"
#include <stdexcept>

//------------------------------------------------------------------------------
// Cheb3D(double t, int N, double Ta, double Tb, Matrix* Cx, Matrix* Cy, Matrix* Cz)
//------------------------------------------------------------------------------
/**
 * Computes the position vector of a celestial body using Chebyshev approximation.
 *
 * This function calculates the position vector of a celestial body at a given time
 * using Chebyshev approximation. It evaluates the Chebyshev polynomial series for
 * each coordinate component (x, y, z) separately and returns the resulting position
 * vector.
 *
 * @param t Time at which to compute the position (in the same units as Ta and Tb).
 * @param N Degree of the Chebyshev polynomial approximation.
 * @param Ta Start time of the approximation interval.
 * @param Tb End time of the approximation interval.
 * @param Cx Pointer to the matrix storing the Chebyshev coefficients for the x-component.
 * @param Cy Pointer to the matrix storing the Chebyshev coefficients for the y-component.
 * @param Cz Pointer to the matrix storing the Chebyshev coefficients for the z-component.
 * @return Position vector of the celestial body at time t.
 */
//------------------------------------------------------------------------------
Matrix Cheb3D(double t, int N, double Ta, double Tb, Matrix* Cx, Matrix* Cy, Matrix* Cz);

#endif // CHEB3D_H
