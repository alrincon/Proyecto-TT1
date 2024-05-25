#include "../include/Cheb3D.h"

//------------------------------------------------------------------------------
// Cheb3D(double t, int N, double Ta, double Tb, Matrix* Cx, Matrix* Cy, Matrix* Cz)
//------------------------------------------------------------------------------
/**
 * Computes the position vector of a celestial body using Chebyshev approximation.
 *
 * This function calculates the position vector of a celestial body at a given time
 * using Chebyshev approximation. It evaluates the Chebyshev polynomial series for
 * each coordinate component (x, y, z) separately and returns the resulting position
 * vector.
 *
 * @param t Time at which to compute the position (in the same units as Ta and Tb).
 * @param N Degree of the Chebyshev polynomial approximation.
 * @param Ta Start time of the approximation interval.
 * @param Tb End time of the approximation interval.
 * @param Cx Pointer to the matrix storing the Chebyshev coefficients for the x-component.
 * @param Cy Pointer to the matrix storing the Chebyshev coefficients for the y-component.
 * @param Cz Pointer to the matrix storing the Chebyshev coefficients for the z-component.
 * @return Position vector of the celestial body at time t.
 */
//------------------------------------------------------------------------------
Matrix Cheb3D(double t, int N, double Ta, double Tb, Matrix* Cx, Matrix* Cy, Matrix* Cz) {
    Matrix ChebApp(1, 3);

    // Check validity
    if ((t < Ta) || (Tb < t)){
        cout << "ERROR: Time out of range in Cheb3D" << endl;
    }

    // Clenshaw algorithm
    double tau = (2*t-Ta-Tb)/(Tb-Ta);

    Matrix f1(1,3);
    Matrix f2(1,3);

    Matrix old_f1(1,3);


    Matrix matAux(1,3);
    for (int i = N; i >= 2; --i) {


        matAux(1,1) = (*Cx)(1,i);
        matAux(1,2) = (*Cy)(1,i);
        matAux(1,3) = (*Cz)(1,i);

        old_f1 = f1;
        f1 = f1*(2.0*tau)-f2+matAux;
        f2 = old_f1;
    }

    matAux(1,1) = (*Cx)(1,1);
    matAux(1,2) = (*Cy)(1,1);
    matAux(1,3) = (*Cz)(1,1);

    ChebApp = f1*tau-f2+matAux;

    return ChebApp;
}
