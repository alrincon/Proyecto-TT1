#ifndef PROYECTO_ANGLESDR_H
#define PROYECTO_ANGLESDR_H

#include "Matrix.h"
#include "SAT_const.h"
#include <stdio.h>
#include <string.h>
#include "Geodetic.h"
#include "LTC.h"
#include "IERS.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "dotProduct.h"
#include "doubler.h"

//------------------------------------------------------------------------------
// void anglesdr(double az1, double az2, double az3, double el1, double el2, double el3,
//               double Mjd1, double Mjd2, double Mjd3, Matrix *rsite1, Matrix *rsite2,
//               Matrix *rsite3, Matrix &r2, Matrix &v2)
//------------------------------------------------------------------------------
/**
 * Determines the position and velocity of a satellite using the angles-only method.
 * This function calculates the position and velocity vectors of a satellite from
 * azimuth and elevation observations at three different times.
 *
 * @param <az1> azimuth angle at the first observation [radians]
 * @param <az2> azimuth angle at the second observation [radians]
 * @param <az3> azimuth angle at the third observation [radians]
 * @param <el1> elevation angle at the first observation [radians]
 * @param <el2> elevation angle at the second observation [radians]
 * @param <el3> elevation angle at the third observation [radians]
 * @param <Mjd1> modified Julian date of the first observation
 * @param <Mjd2> modified Julian date of the second observation
 * @param <Mjd3> modified Julian date of the third observation
 * @param <rsite1> pointer to the site position vector at the first observation [m]
 * @param <rsite2> pointer to the site position vector at the second observation [m]
 * @param <rsite3> pointer to the site position vector at the third observation [m]
 * @param <r2> reference to the resulting position vector of the satellite at the second observation [m]
 * @param <v2> reference to the resulting velocity vector of the satellite at the second observation [m/s]
 * @return none
 * @exception none
 * @note Assumes the input matrices and external data are correctly initialized
 */
//------------------------------------------------------------------------------
void anglesdr (double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, Matrix *rsite1, Matrix *rsite2, Matrix *rsite3, Matrix &r2, Matrix &v2);

#endif
