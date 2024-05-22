#ifndef AZELPA_H
#define AZELPA_H

#include <cmath>
#include "../include/Matrix.h"
#include "../include/dotProduct.h"

void AzElPa(Matrix* s, double& Az, double& El, Matrix& dAds, Matrix& dEds);

#endif // AZELPA_H
