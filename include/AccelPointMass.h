#ifndef ACCELPOINTMASS_H
#define ACCELPOINTMASS_H

#include "../include/Matrix.h"
#include <cmath>
using namespace std;

//------------------------------------------------------------------------------
// Matrix AccelPointMass(Matrix& r, Matrix& s, double GM)
//------------------------------------------------------------------------------
/**
 * Calculates the gravitational acceleration on a satellite due to a point mass.
 * This function computes the acceleration vector of a satellite due to the gravitational
 * attraction of another body, treated as a point mass.
 *
 * @param <r> reference to a Matrix containing the position vector of the satellite
 * @param <s> reference to a Matrix containing the position vector of the point mass
 * @param <GM> gravitational parameter (G * mass) of the point mass
 * @return a Matrix containing the acceleration vector on the satellite
 * @exception none
 * @note no special notes
 */
//------------------------------------------------------------------------------
Matrix AccelPointMass(Matrix& r, Matrix& s, double GM);

//------------------------------------------------------------------------------
// Matrix AccelPointMassT(Matrix& r, Matrix& s, double GM)
//------------------------------------------------------------------------------
/**
 * Calculates the gravitational acceleration on a satellite due to a transposed point mass vector.
 * This function transposes the point mass position vector and then computes the gravitational
 * acceleration on the satellite due to this point mass.
 *
 * @param <r> reference to a Matrix containing the position vector of the satellite
 * @param <s> reference to a Matrix containing the position vector of the point mass
 * @param <GM> gravitational parameter (G * mass) of the point mass
 * @return a Matrix containing the acceleration vector on the satellite
 * @exception none
 * @note no special notes
 */
//------------------------------------------------------------------------------
Matrix AccelPointMassT(Matrix& r, Matrix& s, double GM);

#endif
