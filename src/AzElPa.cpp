#include "../include/AzElPa.h"

//------------------------------------------------------------------------------
// AzElPa(Matrix* s, double& Az, double& El, Matrix& dAds, Matrix& dEds)
//------------------------------------------------------------------------------
/**
 * Calculates azimuth and elevation angles from a given position vector.
 *
 * This function computes the azimuth and elevation angles from a given position
 * vector. It also calculates the partial derivatives of azimuth and elevation
 * with respect to the position vector components.
 *
 * @param s Pointer to the position vector.
 * @param Az Reference to store the azimuth angle (in radians).
 * @param El Reference to store the elevation angle (in radians).
 * @param dAds Matrix to store the partial derivatives of azimuth with respect to
 *             the position vector components.
 * @param dEds Matrix to store the partial derivatives of elevation with respect to
 *             the position vector components.
 */
//------------------------------------------------------------------------------
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

