#include "../include/MeasUpdate.h"

//Los vectores son (n,1)
void MeasUpdate(Matrix* z, Matrix* g, Matrix* s, Matrix*G , int n, Matrix& K, Matrix& x, Matrix& P){
    int m = z->getFilas();
    Matrix Inv_W(m,m);

    for(int i = 1; i <= m; i++) {
        Inv_W(i, i) = (*s)(i, 1) * (*s)(i ,1);    // Inverse weight(measurement covariance)
    }

    // Kalman gain
    Matrix prod = P*(*G).transpose()*(Inv_W+(*G)*P*(*G).transpose()).inverse();
    (K).redefine(&prod);

    // State update

    x = x + ((*z)-(*g))*K.transpose();

    Matrix id(n,n, true);
    // Covariance update
    P = (id-K*(*G))*P;

}
//[K, x, P] = MeasUpdate(x, z, g, s, G, P, n)
void MeasUpdate(double  z, double g, double s, Matrix*G , int n, Matrix& K, Matrix& x, Matrix& P){
    Matrix zM(1,1);
    zM(1,1) = z;

    Matrix gM(1,1);
    gM(1,1) = g;

    Matrix sM(1,1);
    sM(1,1) = s;

    MeasUpdate(&zM, &gM, &sM, G ,n , K, x, P);
}