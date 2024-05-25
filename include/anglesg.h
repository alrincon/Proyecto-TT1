#ifndef PROYECTO_ANGLESG_H
#define PROYECTO_ANGLESG_H

#include "Matrix.h"
#include "Geodetic.h"
#include "LTC.h"
#include "IERS.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "GHAMatrix.h"
#include "PoleMatrix.h"
#include "SAT_const.h"
#include "dotProduct.h"
#include "gibbs.h"
#include "hgibbs.h"
#include "elements.h"
#include "angl.h"
#include "largestRoot.h"
#include "SAT_const.h"

//------------------------------------------------------------------------------
// anglesg(double az1, double az2, double az3, double el1, double el2, double el3,
//         double Mjd1, double Mjd2, double Mjd3, Matrix *Rs1, Matrix *Rs2, Matrix *Rs3,
//         Matrix &r2, Matrix &v2)
//------------------------------------------------------------------------------
/**
 * Calculates the angles between two vectors from two different observers
 * to a common target, using Gauss's method.
 *
 * This function computes the angles between two vectors from two different
 * observers to a common target, based on the observations' azimuth, elevation,
 * and observation times. It uses Gauss's method to estimate the angles and
 * determines the position and velocity of the target at the second observation
 * time.
 *
 * @param az1 Azimuth of the target as observed from the first observer (in radians).
 * @param az2 Azimuth of the target as observed from the second observer (in radians).
 * @param az3 Azimuth of the target as observed from the third observer (in radians).
 * @param el1 Elevation of the target as observed from the first observer (in radians).
 * @param el2 Elevation of the target as observed from the second observer (in radians).
 * @param el3 Elevation of the target as observed from the third observer (in radians).
 * @param Mjd1 Modified Julian Date (MJD) of the first observation.
 * @param Mjd2 Modified Julian Date (MJD) of the second observation.
 * @param Mjd3 Modified Julian Date (MJD) of the third observation.
 * @param Rs1 Pointer to the position vector of the first observer (Earth-fixed system).
 * @param Rs2 Pointer to the position vector of the second observer (Earth-fixed system).
 * @param Rs3 Pointer to the position vector of the third observer (Earth-fixed system).
 * @param r2 Position vector of the target at the second observation time (output).
 * @param v2 Velocity vector of the target at the second observation time (output).
 */
//------------------------------------------------------------------------------
void anglesg (double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, Matrix *Rs1, Matrix *Rs2, Matrix *Rs3, Matrix &r2, Matrix &v2);


#endif
