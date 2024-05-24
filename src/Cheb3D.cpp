#include "../include/Cheb3D.h"

Matrix Cheb3D(double t, int N, double Ta, double Tb, Matrix* Cx, Matrix* Cy, Matrix* Cz) {
    Matrix ChebApp(1, 3);
    /*
    // Check validity
    if ((t < Ta) || (Tb < t)) {
        throw runtime_error("ERROR: Time out of range in Cheb3D::Value\n");
    }

    // Clenshaw algorithm
    double tau = (2 * t - Ta - Tb) / (Tb - Ta);


    Matrix f1(1, 3);
    Matrix f2(1, 3);
    double v[3] = {(*Cx)(1,1),(*Cy)(1,1),(*Cz)(1,1)};
    Matrix aux(1,3,v,3);

    for (int i = N; i >= 2; i--) {
        Matrix old_f1(1,3);
        old_f1 = f1;
        double w[3] = {(*Cx)(1,i),(*Cy)(1,i),(*Cz)(1,i)};
        Matrix aux2(1,3,w,3);
        f1 = f1 * tau * 2 - f2 + aux2;
        f2 = old_f1;
    }
    res = f1 * tau - f2 + aux;

    return res;*/

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
        f1 = f1*(2*tau)-f2+matAux;
        f2 = old_f1;
    }

    matAux(1,1) = (*Cx)(1,1);
    matAux(1,2) = (*Cy)(1,1);
    matAux(1,3) = (*Cz)(1,1);

    ChebApp = f1*tau-f2+matAux;

    return ChebApp;
}
