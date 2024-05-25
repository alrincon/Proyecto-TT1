#ifndef PROYECTO_ANGL_H
#define PROYECTO_ANGL_H

#include "../include/Matrix.h"
#include "../include/sign_.h"
#include <cmath>

//------------------------------------------------------------------------------
// double sign(double x)
//------------------------------------------------------------------------------
/*
 * Returns the sign of a number.
 * This function returns -1 if the input is negative, otherwise it returns 1.
 *
 * @param <x> the input value
 * @return -1 if the input is negative, otherwise 1
 * @exception none
 * @note useful for ensuring the correct sign in calculations
 */
//------------------------------------------------------------------------------
double sign (double x);

//------------------------------------------------------------------------------
// double angl(Matrix* vec1, Matrix* vec2)
//------------------------------------------------------------------------------
/**
 * Calculates the angle between two vectors.
 * This function computes the angle between two vectors in radians using the dot product.
 *
 * @param <vec1> pointer to the first vector (Matrix)
 * @param <vec2> pointer to the second vector (Matrix)
 * @return the angle between the two vectors in radians, or a large value if undefined
 * @exception none
 * @note returns a large value if either vector is zero-length
 */
//------------------------------------------------------------------------------
double angl (Matrix* vec1, Matrix* vec2);

#endif
