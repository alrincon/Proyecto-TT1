#include "../include/AccelPointMass.h"
#include <cmath>

/*
%--------------------------------------------------------------------------
%
% AccelPointMass: Computes the perturbational acceleration due to a point
%				  mass
%
% Inputs:
%   r           Satellite position vector
%   s           Point mass position vector
%   GM          Gravitational coefficient of point mass
%
% Output:
%   a    		Acceleration (a=d^2r/dt^2)
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
*/

double AccelPointMass(double r, double s, double GM){
    //Relative position vector of satellite w.r.t. point mass
    double d = r - s;

    //Acceleration
    return -GM * ( d/(abs(d)^3) + s/(abs(s)^3) );
}