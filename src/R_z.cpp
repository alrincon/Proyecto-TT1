#include "../include/Matrix.h"
#include "../include/R_z.h"
#include <cmath>

//------------------------------------------------------------------------------
// Matrix R_z(double angle)
//------------------------------------------------------------------------------
/**
 * Computes the rotation matrix about the z-axis.
 *
 * @param angle  Angle of rotation (radians).
 * @return       Rotation matrix about the z-axis.
 */
//------------------------------------------------------------------------------
Matrix R_z(double angle){
    Matrix rotmat(3,3);

    double C = cos(angle);
    double S = sin(angle);

    rotmat(1,1) =      C;  rotmat(1,2) =   S;  rotmat(1,3) = 0.0;
    rotmat(2,1) = -1.0*S;  rotmat(2,2) =   C;  rotmat(2,3) = 0.0;
    rotmat(3,1) =    0.0;  rotmat(3,2) = 0.0;  rotmat(3,3) = 1.0;

    return rotmat;
}

