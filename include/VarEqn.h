#ifndef PROYECTO_VAREQN_H
#define PROYECTO_VAREQN_H

#include "Matrix.h"
#include "SAT_const.h"
#include "IERS.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "PoleMatrix.h"
#include "NutMatrix.h"
#include "GHAMatrix.h"
#include "AccelHarmonic.h"
#include "G_AccelHarmonic.h"
#include "types.h"

/**
 * @brief Computes the time derivative of the combined state vector and state transition matrix.
 *
 * This function calculates the time derivative of the combined state vector and state transition matrix
 * for the given time 'x' and the provided state vector 'yPhi'. It utilizes various auxiliary parameters
 * and functions to perform the computations.
 *
 * @param x Time at which to compute the time derivative, given in seconds.
 * @param yPhi Pointer to the Matrix object representing the combined state vector and state transition matrix.
 *             The state vector should be arranged as [r, v, Phi], where 'r' is the position vector,
 *             'v' is the velocity vector, and 'Phi' is the state transition matrix.
 * @return Matrix A Matrix object containing the time derivative of the combined state vector and state transition matrix.
 *
 * @note The function internally uses various auxiliary parameters and functions, such as IERS(), timediff(),
 *       PrecMatrix(), NutMatrix(), PoleMatrix(), GHAMatrix(), AccelHarmonic(), and G_AccelHarmonic(), to compute
 *       the required values and matrices.
 */
Matrix VarEqn(double x, Matrix* yPhi);

Matrix VarEqnOUT(double x, Matrix& yPhi);

#endif
