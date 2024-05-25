#ifndef PROYECTO_TIMEUPDATE_H
#define PROYECTO_TIMEUPDATE_H

#include "Matrix.h"

/**
 * @brief Updates the covariance matrix of the system's state estimate incorporating process noise.
 *
 * This function updates the covariance matrix (P) of the state estimate using the state
 * transition matrix (Phi) and the process noise covariance matrix (Qdt). It is typically used
 * in the predict step of Kalman filtering, where the uncertainty in the state estimate
 * increases due to the process dynamics and external influences modeled by the process noise.
 *
 * @param P The covariance matrix of the state estimate that will be updated.
 * @param Phi Pointer to the state transition matrix, which defines how the state evolves from
 *            one time step to the next without considering the process noise.
 * @param Qdt Pointer to the process noise covariance matrix, which adds uncertainty from
 *            the model and environment.
 */
void TimeUpdate(Matrix& P, Matrix* Phi, Matrix* Qdt);

/**
 * @brief Updates the covariance matrix of the system's state estimate without incorporating process noise.
 *
 * This overload of the TimeUpdate function updates the covariance matrix (P) using only the state
 * transition matrix (Phi). It is suitable for cases where the model's process noise is negligible
 * or is being modeled separately. This function focuses on the uncertainty propagation through
 * the system dynamics alone.
 *
 * @param P The covariance matrix of the state estimate that will be updated.
 * @param Phi Pointer to the state transition matrix, which defines how the state evolves from
 *            one time step to the next.
 */
void TimeUpdate(Matrix& P, Matrix* Phi);

#endif
