#include "../include/Matrix.h"
#include "../include/R_x.h"
#include <cmath>

//------------------------------------------------------------------------------
// Matrix R_x(double angle)
//------------------------------------------------------------------------------
/**
 * Computes the rotation matrix about the x-axis.
 *
 * @param angle  Angle of rotation (radians).
 * @return       Rotation matrix about the x-axis.
 */
//------------------------------------------------------------------------------
Matrix R_x(double angle){
    Matrix rotmat(3,3);

    double C = cos(angle);
    double S = sin(angle);

    rotmat(1,1) = 1.0;  rotmat(1,2) =    0.0;  rotmat(1,3) = 0.0;
    rotmat(2,1) = 0.0;  rotmat(2,2) =      C;  rotmat(2,3) =   S;
    rotmat(3,1) = 0.0;  rotmat(3,2) = -1.0*S;  rotmat(3,3) =   C;

    return rotmat;
}

