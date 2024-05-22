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

//los vectores son n,1
void AzElPa(Matrix* s, double& Az, double& El, Matrix& dAds, Matrix& dEds) {
    double pi2 = 2.0 * M_PI;

    double rho = sqrt((*s)(1,1)*(*s)(1,1)+(*s)(2,1)*(*s)(2,1));

    // Angles
    Az = atan2((*s)(1,1),(*s)(2,1));

    if (Az<0.0) {
        Az = Az + pi2;
    }

    El = atan ( (*s)(3,1) / rho );

    // Partials
    dAds(1,1) = (*s)(2,1)/(rho*rho);
    dAds(2,1) = -(*s)(1,1)/(rho*rho);
    dAds(3,1) = 0;

    dEds(1,1) = -(*s)(1,1)*(*s)(3,1)/(rho*dotProduct(s,s));
    dEds(2,1) = -(*s)(2,1)*(*s)(3,1)/(rho*dotProduct(s,s));
    dEds(3,1) = rho/dotProduct(s,s);
}

