#include "../include/AccelPointMass.h"
#include <cmath>
#include "../include/Vector.h"

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

double AccelPointMass(Vector r, Vector s, double GM){
    //Relative position vector of satellite w.r.t. point mass
    Vector d = r - s;

    //Acceleration
    return -GM * ( d/(d.norm()) + s/(s.norm()^3) );
}