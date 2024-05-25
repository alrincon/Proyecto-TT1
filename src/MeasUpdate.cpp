#include "../include/MeasUpdate.h"

//------------------------------------------------------------------------------
// void MeasUpdate(Matrix* z, Matrix* g, Matrix* s, Matrix* G, int n, Matrix& K, Matrix& x, Matrix& P)
//------------------------------------------------------------------------------
/**
 * Updates the state estimate and covariance matrix based on a measurement.
 *
 * @param z Measurement vector.
 * @param g Predicted measurement vector.
 * @param s Measurement noise vector.
 * @param G Measurement matrix.
 * @param n Dimension of the state vector.
 * @param K Kalman gain matrix (output).
 * @param x State vector to be updated (input/output).
 * @param P Covariance matrix to be updated (input/output).
 */
//------------------------------------------------------------------------------
void MeasUpdate(Matrix* z, Matrix* g, Matrix* s, Matrix*G , int n, Matrix& K, Matrix& x, Matrix& P){
    int m = z->getFilas();
    Matrix Inv_W(m,m);

    Matrix xT = x.transpose();

    for(int i = 1; i <= m; i++) {
        Inv_W(i, i) = (*s)(i, 1) * (*s)(i ,1);    // Inverse weight(measurement covariance)
    }

    // Kalman gain
    Matrix prod = P*(*G).transpose()*(Inv_W+(*G)*P*(*G).transpose()).inverse();
    (K).redefine(&prod);

    // State update

    xT = xT + K*((*z)-(*g));

    x = xT.transpose();

    Matrix id(n,n, true);
    // Covariance update
    P = (id-K*(*G))*P;

}

//------------------------------------------------------------------------------
// void MeasUpdate(double z, double g, double s, Matrix* G, int n, Matrix& K, Matrix& x, Matrix& P)
//------------------------------------------------------------------------------
/**
 * Updates the state estimate and covariance matrix based on a single scalar measurement.
 * This function is a simplified version of MeasUpdate() for scalar measurements.
 *
 * @param z Measurement scalar.
 * @param g Predicted measurement scalar.
 * @param s Measurement noise scalar.
 * @param G Measurement matrix.
 * @param n Dimension of the state vector.
 * @param K Kalman gain matrix (output).
 * @param x State vector to be updated (input/output).
 * @param P Covariance matrix to be updated (input/output).
 */
//------------------------------------------------------------------------------
void MeasUpdate(double  z, double g, double s, Matrix*G , int n, Matrix& K, Matrix& x, Matrix& P){
    Matrix zM(1,1);
    zM(1,1) = z;

    Matrix gM(1,1);
    gM(1,1) = g;

    Matrix sM(1,1);
    sM(1,1) = s;

    MeasUpdate(&zM, &gM, &sM, G ,n , K, x, P);
}