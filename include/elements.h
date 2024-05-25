#ifndef PROYECTO_ELEMENTS_H
#define PROYECTO_ELEMENTS_H

#include "SAT_const.h"
#include "Matrix.h"
#include "crossProduct.h"
#include "dotProduct.h"
#include <cmath>
#include "realmod.h"

//------------------------------------------------------------------------------
// elements(Matrix* y, double &p, double &a, double &e, double &i, double &Omega, double &omega, double &M)
//------------------------------------------------------------------------------
/**
 * Computes orbital elements from state vector.
 *
 * This function computes orbital elements (semi-latus rectum, semi-major axis,
 * eccentricity, inclination, longitude of the ascending node, argument of
 * perigee, and mean anomaly) from the state vector (position and velocity).
 *
 * @param y Pointer to the state vector (position and velocity).
 * @param p Semi-latus rectum.
 * @param a Semi-major axis.
 * @param e Eccentricity.
 * @param i Inclination.
 * @param Omega Longitude of the ascending node.
 * @param omega Argument of perigee.
 * @param M Mean anomaly.
 */
//------------------------------------------------------------------------------
void elements (Matrix* y, double &p, double &a, double &e, double &i, double &Omega, double &omega, double &M);

#endif
