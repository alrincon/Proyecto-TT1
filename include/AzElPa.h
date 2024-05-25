#ifndef AZELPA_H
#define AZELPA_H

#include <cmath>
#include "../include/Matrix.h"
#include "../include/dotProduct.h"

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
void AzElPa(Matrix* s, double& Az, double& El, Matrix& dAds, Matrix& dEds);

#endif // AZELPA_H
