#include "../include/AccelPointMass.h"

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
Matrix AccelPointMass(Matrix& r, Matrix& s, double GM){
    //Relative position vector of satellite w.r.t. point mass
    Matrix d = r - s;

    //Acceleration
    return (d*(1.0/(pow(d.norm(),3))) + s*(1.0/(pow(s.norm(),3))))*(-GM);
}

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
Matrix AccelPointMassT(Matrix& r, Matrix& s, double GM){
    Matrix sT(r.getFilas(),1);

    for(int i = 1; i < r.getFilas(); i++){
        sT(i,1) = s(1,i);
    }

    return AccelPointMass(r, sT, GM);
}