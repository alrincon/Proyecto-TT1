#include "../include/MeasUpdate.h"

void MeasUpdate(Matrix* z, Matrix* g, Matrix* s, Matrix*G , int n, Matrix& K, Matrix& x, Matrix& P){
    int m = z->getColumnas();
    Matrix Inv_W(m,m);

    for(int i = 1; i <= m; i++) {
        Inv_W(i, i) = (*s)(1, i) * (*s)(1 ,i);    // Inverse weight(measurement covariance)
    }

    // Kalman gain
    K = P*(*G).transpose()*(Inv_W+(*G)*P*(*G).transpose()).inverse();

    // State update
    x = x + K*(z-g);

    Matrix id(n,n, true);
    // Covariance update
    P = (id-K*(*G))*P;

}