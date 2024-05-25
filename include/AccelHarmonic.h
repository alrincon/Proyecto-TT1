#ifndef PROYECTO_ACCELHARMONIC_H
#define PROYECTO_ACCELHARMONIC_H

#include "Matrix.h"
#include "SAT_const.h"
#include "Legendre.h"
#include <cmath>

//------------------------------------------------------------------------------
// Matrix AccelHarmonic(Matrix* r, Matrix* E, int n_max, int m_max)
//------------------------------------------------------------------------------
/**
 * Calculates the acceleration due to the harmonic gravity field.
 * This function computes the gravitational acceleration based on the spherical harmonic
 * expansion of Earth's gravity field.
 *
 * @param <r> pointer to a Matrix containing the position vector in inertial coordinates
 * @param <E> pointer to a Matrix representing the transformation from inertial to body-fixed coordinates
 * @param <n_max> maximum degree of the spherical harmonics
 * @param <m_max> maximum order of the spherical harmonics
 * @return a Matrix containing the acceleration vector in inertial coordinates
 * @exception none
 * @note caller is responsible for managing the memory of the input matrices r and E
 */
//------------------------------------------------------------------------------
Matrix AccelHarmonic(Matrix* r, Matrix* E, int n_max, int m_max);


#endif
