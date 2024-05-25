#include "../include/EccAnom.h"

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
double EccAnom(double M, double e) {
    int maxit = 15;
    int i = 1;

    // Starting value
    M = fmod(M, 2.0 * M_PI);

    double E;
    if (e < 0.8) {
        E = M;
    } else {
        E = M_PI;
    }

    double f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));

    // Iteration
    while (abs(f) > 1e2 *  2.2204e-16) {
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        i++;
        if (i == maxit) {
            throw runtime_error("convergence problems in EccAnom");
        }
    }

    return E;
}

