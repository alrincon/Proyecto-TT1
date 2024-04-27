#include "../include/Cheb3D.h"

Matrix Cheb3D(double t, int N, double Ta, double Tb, Matrix& Cx, Matrix& Cy, Matrix& Cz) {
    // Check validity
    if ((t < Ta) || (Tb < t)) {
        throw runtime_error("ERROR: Time out of range in Cheb3D::Value\n");
    }

    // Clenshaw algorithm
    double tau = (2 * t - Ta - Tb) / (Tb - Ta);

    Matrix f1(1, 3);
    Matrix f2(1, 3);
    double v[3] = {Cx(1,1),Cy(1,1),Cz(1,1)};
    Matrix aux(1,3,v,3);

    for (int i = N - 1; i >= 1; --i) {
        Matrix old_f1 = f1;
        double w[3] = {Cx(i,1),Cy(i,1),Cz(i,1)};
        Matrix aux2(1,3,w,3);
        f1 = f1 * tau * 2 - f2 + aux2;

        f2 = old_f1;
    }

    return f1 * tau - f2 + aux;
}
