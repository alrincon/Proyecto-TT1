#include "../include/Matrix.h"
#include "../include/R_x.h"
#include <cmath>

/*

%--------------------------------------------------------------------------
%  input:
%    angle       - angle of rotation [rad]
%
%  output:
%    rotmat      - vector result
%--------------------------------------------------------------------------
function [rotmat] = R_z(angle)

C = cos(angle);
S = sin(angle);
rotmat = zeros(3,3);

rotmat(1,1) =      C;  rotmat(1,2) =   S;  rotmat(1,3) = 0.0;
rotmat(2,1) = -1.0*S;  rotmat(2,2) =   C;  rotmat(2,3) = 0.0;
rotmat(3,1) =    0.0;  rotmat(3,2) = 0.0;  rotmat(3,3) = 1.0;


 */

Matrix R_x(double angle){
    Matrix sol(3,3);

    double C = cos(angle);
    double S = sin(angle);

    sol(1,1) = 1.0;
    sol(1,2) = 0.0;
    sol(1,3) = 0.0;

    sol(2,1) = 0.0;
    sol(2,2) = C;
    sol(2,3) = S;

    sol(3,1) = 0.0;
    sol(3,2) = -1.0*S;
    sol(3,3) = C;

    return sol;
}

