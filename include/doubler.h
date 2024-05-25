#ifndef PROYECTO_DOUBLER_H
#define PROYECTO_DOUBLER_H

#include <cmath>
#include "SAT_const.h"
#include "Matrix.h"
#include "dotProduct.h"
#include "crossProduct.h"
#include <stdio.h>
#include <string.h>

//------------------------------------------------------------------------------
// doubler(double cc1, double cc2, double magrsite1, double magrsite2, double magr1in, double magr2in, Matrix *los1, Matrix *los2, Matrix *los3, Matrix *rsite1, Matrix *rsite2, Matrix *rsite3, double t1, double t3, char direct, Matrix &r2, Matrix &r3, double &f1, double &f2, double &q1, double &magr1, double &magr2, double &a, double &deltae32)
//------------------------------------------------------------------------------
/**
 * Computes the vectors r2 and r3, as well as other parameters, for a double r approach.
 *
 * This function computes the vectors r2 and r3 for a double r approach, as well as other
 * parameters such as f1, f2, q1, magr1, magr2, a, and deltae32. It accepts several input
 * parameters and modifies some of them by reference.
 *
 * @param cc1 Double r coefficient for the first site.
 * @param cc2 Double r coefficient for the second site.
 * @param magrsite1 Magnitude of the position vector of the first site.
 * @param magrsite2 Magnitude of the position vector of the second site.
 * @param magr1in Magnitude of the input position vector r1.
 * @param magr2in Magnitude of the input position vector r2.
 * @param los1 Pointer to the line of sight unit vector for the first site.
 * @param los2 Pointer to the line of sight unit vector for the second site.
 * @param los3 Pointer to the line of sight unit vector for the third site.
 * @param rsite1 Pointer to the position vector of the first site.
 * @param rsite2 Pointer to the position vector of the second site.
 * @param rsite3 Pointer to the position vector of the third site.
 * @param t1 Time at the first site.
 * @param t3 Time at the third site.
 * @param direct Direct or indirect double r approach ('y' or 'n').
 * @param r2 Computed position vector r2.
 * @param r3 Computed position vector r3.
 * @param f1 Resulting time difference f1.
 * @param f2 Resulting time difference f2.
 * @param q1 Resulting time difference q1.
 * @param magr1 Magnitude of the computed position vector r1.
 * @param magr2 Magnitude of the computed position vector r2.
 * @param a Semi-major axis of the computed orbit.
 * @param deltae32 Difference between eccentric anomaly of orbits 3 and 2.
 *
 * @note The function modifies the values of r2, r3, f1, f2, q1, magr1, magr2, a, and deltae32.
 */
//------------------------------------------------------------------------------
void doubler(double cc1, double cc2, double magrsite1, double magrsite2, double magr1in, double magr2in, Matrix *los1, Matrix *los2, Matrix *los3, Matrix *rsite1, Matrix *rsite2, Matrix *rsite3, double t1, double t3, char direct, Matrix &r2, Matrix &r3, double &f1, double &f2, double &q1, double &magr1, double &magr2, double &a, double &deltae32);

#endif
