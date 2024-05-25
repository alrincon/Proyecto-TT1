#ifndef PROYECTO_G_ACCELHARMONIC_H
#define PROYECTO_G_ACCELHARMONIC_H

#include "Matrix.h"
#include "AccelHarmonic.h"

//------------------------------------------------------------------------------
// G_AccelHarmonic(Matrix* r, Matrix* U, int n_max, int m_max)
//------------------------------------------------------------------------------
/**
 * Computes the gradient of the gravitational acceleration due to harmonic coefficients.
 *
 * This function calculates the gradient of the gravitational acceleration due
 * to the harmonic coefficients at a given position.
 *
 * @param r Pointer to the position vector.
 * @param U Pointer to the harmonic coefficients.
 * @param n_max Maximum degree of harmonic coefficients.
 * @param m_max Maximum order of harmonic coefficients.
 * @return Gradient of the gravitational acceleration.
 */
//------------------------------------------------------------------------------
Matrix G_AccelHarmonic(Matrix* r, Matrix* U, int n_max, int m_max);

#endif
