#ifndef PROYECTO_GIBBS_H
#define PROYECTO_GIBBS_H

#include "Matrix.h"
#include "SAT_const.h"
#include "unit.h"
#include "crossProduct.h"
#include "dotProduct.h"
#include "angl.h"
#include <stdio.h>
#include <string.h>

//------------------------------------------------------------------------------
// gibbs(Matrix *r1, Matrix *r2, Matrix *r3, Matrix &v2, double &theta,
//       double &theta1, double &copa, char* &error)
//------------------------------------------------------------------------------
/**
 * Computes the velocity vector at the second point using the Gibbs method
 * given the position vectors at three points in space.
 *
 * @param r1 Pointer to the position vector at the first point.
 * @param r2 Pointer to the position vector at the second point.
 * @param r3 Pointer to the position vector at the third point.
 * @param v2 Reference to the velocity vector at the second point (output).
 * @param theta Reference to the first angular change (output).
 * @param theta1 Reference to the second angular change (output).
 * @param copa Reference to the angle between the plane of the three points
 *             and the first position vector (output).
 * @param error Reference to the error message (output).
 */
//------------------------------------------------------------------------------
void gibbs(Matrix *r1, Matrix *r2, Matrix *r3, Matrix &v2, double &theta, double &theta1, double &copa, char* &error);

#endif
