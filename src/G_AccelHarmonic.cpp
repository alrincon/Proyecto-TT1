#include "../include/G_AccelHarmonic.h"

Matrix G_AccelHarmonic(Matrix* r, Matrix* U, double n_max, double m_max){
    double d = 1.0;   // Position increment [m]

    Matrix G(3,3);
    Matrix dr(3,1);

    // Gradient
    for(int i = 1; i <= 3; i++){
            // Set offset in i-th component of the position vector
        dr(i,1) = d;
        // Acceleration difference
        Matrix a = *r+dr*(1/2);
        Matrix b = *r-dr*(1/2);
        Matrix da = AccelHarmonic ( &a,U, n_max, m_max ) - AccelHarmonic ( &b,U, n_max, m_max );

        // Derivative with respect to i-th axis
        for(int j = 1; j <= 3; j++){
            G(j,i) = da(j,1)/d;
        }
    }

    return G;
}