#ifndef PROYECTO_HGIBBS_H
#define PROYECTO_HGIBBS_H

#include "Matrix.h"
#include "SAT_const.h"
#include <stdio.h>
#include <string.h>
#include "crossProduct.h"
#include "dotProduct.h"
#include "unit.h"

//------------------------------------------------------------------------------
// void hgibbs (Matrix* r1, Matrix* r2, Matrix* r3,
//              double Mjd1, double Mjd2, double Mjd3,
//              Matrix& v2, double &theta, double &theta1,
//              double &copa, char* &error)
//------------------------------------------------------------------------------
/**
 * Implements the Herrick-Gibbs method for solving the three-body problem,
 * providing the velocity vector at the second position based on three
 * position vectors and their respective times.
 *
 * @param r1    Pointer to the position vector at the first time.
 * @param r2    Pointer to the position vector at the second time.
 * @param r3    Pointer to the position vector at the third time.
 * @param Mjd1  Modified Julian Date (MJD) of the first time.
 * @param Mjd2  Modified Julian Date (MJD) of the second time.
 * @param Mjd3  Modified Julian Date (MJD) of the third time.
 * @param v2    Output velocity vector at the second position.
 * @param theta Angle between r1 and r2.
 * @param theta1 Angle between r2 and r3.
 * @param copa  Angle between the normal to the plane of r2 and r3 and r1.
 * @param error Error message if any.
 */
//------------------------------------------------------------------------------
void hgibbs (Matrix* r1, Matrix* r2, Matrix* r3, double Mjd1, double Mjd2, double Mjd3, Matrix& v2, double &theta, double &theta1,double &copa, char* &error);

#endif
