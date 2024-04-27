#include "../include/EccAnom.h"

/*%--------------------------------------------------------------------------
%
% Purpose:
%   Computes the eccentric anomaly for elliptic orbits
%
% Inputs:
%   M         Mean anomaly in [rad]
%   e         Eccentricity of the orbit [0,1]
%
% Output:
%             Eccentric anomaly in [rad]
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------*/

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
    while (abs(f) > 1e2 * numeric_limits<double>::epsilon()) {
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        i++;
        if (i == maxit) {
            throw runtime_error("convergence problems in EccAnom");
        }
    }

    return E;
}

