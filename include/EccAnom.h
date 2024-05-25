#ifndef ECCANOM_H
#define ECCANOM_H

#include <cmath>
#include <stdexcept>
#include <limits>

using namespace std;

//------------------------------------------------------------------------------
// EccAnom(double M, double e)
//------------------------------------------------------------------------------
/**
 * Computes the eccentric anomaly for a given mean anomaly and eccentricity.
 *
 * This function computes the eccentric anomaly for a given mean anomaly and
 * eccentricity using the Newton-Raphson method. It iterates until convergence
 * or until the maximum number of iterations is reached.
 *
 * @param M Mean anomaly.
 * @param e Eccentricity.
 *
 * @return Eccentric anomaly.
 *
 * @throws runtime_error if convergence problems occur during iteration.
 */
//------------------------------------------------------------------------------
double EccAnom(double M, double e);

#endif // ECCANOM_H
