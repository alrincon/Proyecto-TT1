#ifndef PROYECTO_MEASUPDATE_H
#define PROYECTO_MEASUPDATE_H

#include "Matrix.h"

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
void MeasUpdate(Matrix* z, Matrix* g, Matrix* s, Matrix*G , int n, Matrix& K, Matrix& x, Matrix& P);

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
void MeasUpdate(double z, double g, double s, Matrix*G , int n, Matrix& K, Matrix& x, Matrix& P);

#endif
