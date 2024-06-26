#include "../include/Matrix.h"
#include "../include/R_y.h"
#include <cmath>

//------------------------------------------------------------------------------
// Matrix R_y(double angle)
//------------------------------------------------------------------------------
/**
 * Computes the rotation matrix about the y-axis.
 *
 * @param angle  Angle of rotation (radians).
 * @return       Rotation matrix about the y-axis.
 */
//------------------------------------------------------------------------------
Matrix R_y(double angle){
    Matrix rotmat(3,3);

    double C = cos(angle);
    double S = sin(angle);

    rotmat(1,1) =   C;  rotmat(1,2) = 0.0;  rotmat(1,3) = -1.0*S;
    rotmat(2,1) = 0.0;  rotmat(2,2) = 1.0;  rotmat(2,3) =    0.0;
    rotmat(3,1) =   S;  rotmat(3,2) = 0.0;  rotmat(3,3) =      C;

    return rotmat;
}