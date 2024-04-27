#include "../include/AzElPa.h"

//--------------------------------------------------------------------------
//
// Purpose:
//  Computes azimuth, elevation and partials from local tangent coordinates
//
// Input:
//   s      Topocentric local tangent coordinates (East-North-Zenith frame)
//
// Outputs:
//   A      Azimuth [rad]
//   E      Elevation [rad]
//   dAds   Partials of azimuth w.r.t. s
//   dEds   Partials of elevation w.r.t. s
//
// Last modified:   2015/08/12   M. Mahooti
//
//--------------------------------------------------------------------------


void AzElPa(const double s[3], double& Az, double& El, double dAds[3], double dEds[3]) {
    const double pi2 = 2.0 * M_PI;

    double rho = std::sqrt(s[0] * s[0] + s[1] * s[1]);

    // Angles
    Az = std::atan2(s[0], s[1]);
    if (Az < 0.0) {
        Az += pi2;
    }

    El = std::atan(s[2] / rho);

    // Partials
    dAds[0] = s[1] / (rho * rho);
    dAds[1] = -s[0] / (rho * rho);
    dAds[2] = 0.0;

    double dot_s = s[0] * s[0] + s[1] * s[1] + s[2] * s[2];
    dEds[0] = -s[0] * s[2] / rho;
    dEds[1] = -s[1] * s[2] / rho;
    dEds[2] = rho / dot_s;
}

