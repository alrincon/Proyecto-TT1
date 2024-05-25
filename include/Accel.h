#ifndef PROYECTO_ACCEL_H
#define PROYECTO_ACCEL_H

#include "Matrix.h"
#include "IERS.h"
#include "timediff.h"
#include "SAT_const.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "Mjday_TDB.h"
#include "AccelHarmonic.h"
#include "AccelPointMass.h"
#include "JPL_Eph_DE430.h"
#include "types.h"

//------------------------------------------------------------------------------
// Matrix Accel(double x, Matrix* Y)
//------------------------------------------------------------------------------
/**
 * Calculates the acceleration of a spacecraft based on various perturbations.
 * This function computes the acceleration due to harmonic gravity field, luni-solar
 * perturbations, and planetary perturbations.
 *
 * @param <x> time since epoch in seconds
 * @param <Y> pointer to a Matrix containing state vector [position; velocity]
 * @return a Matrix containing the derivative of the state vector [velocity; acceleration]
 * @exception none
 * @note caller is responsible for managing the memory of the input Matrix Y
 */
//------------------------------------------------------------------------------
Matrix Accel(double x, Matrix* Y);

//------------------------------------------------------------------------------
// Matrix AccelOUT(double x, Matrix& Y)
//------------------------------------------------------------------------------
/**
 * Wrapper function for Accel.
 * This function provides a more convenient interface for calling Accel with a Matrix reference.
 *
 * @param <x> time since epoch in seconds
 * @param <Y> reference to a Matrix containing state vector [position; velocity]
 * @return a Matrix containing the derivative of the state vector [velocity; acceleration]
 * @exception none
 * @note no special notes
 */
//------------------------------------------------------------------------------
Matrix AccelOUT(double x, Matrix& Y);
#endif
