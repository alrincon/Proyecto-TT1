#include "../include/Geodetic.h"
#include "../include/SAT_const.h"

//------------------------------------------------------------------------------
// Geodetic(double& lon, double& lat, double& h, Matrix* r)
//------------------------------------------------------------------------------
/**
 * Computes the geodetic coordinates (longitude, latitude, and altitude) from
 * the Cartesian coordinates (X, Y, Z).
 *
 * @param lon Longitude (output).
 * @param lat Latitude (output).
 * @param h Altitude (output).
 * @param r Cartesian coordinates (input).
 */
//------------------------------------------------------------------------------
void Geodetic(double& lon, double& lat, double& h, Matrix* r){

    double eps = 0.0000001;
    double R_equ = R_Earth;
    double f     = f_Earth;

    double epsRequ = eps*R_equ;        //Convergence criterion
    double e2      = f*(2.0-f);        // Square of eccentricity

    double X = (*r)(1,1);                   // Cartesian coordinates
    double Y = (*r)(2,1);
    double Z = (*r)(3,1);
    double rho2 = X*X + Y*Y;           // Square of distance from z-axis


                                                            // Check validity of input data
    if ((*r).norm() < 10e-6) {
        lon = 0.0;
        lat = 0.0;
        h = -R_Earth;
    }

    // Iteration
    double dZ = e2*Z;
    double ZdZ;
    double Nh;
    double N;

    while(1) {
        ZdZ = Z + dZ;
        Nh = sqrt(rho2 + ZdZ * ZdZ);
        double SinPhi = ZdZ / Nh;                    // Sine of geodetic latitude
        N = R_equ / sqrt(1.0 - e2 * SinPhi * SinPhi);
        double dZ_new = N * e2 * SinPhi;
        if (fabs(dZ - dZ_new) < epsRequ) {
            break;
        }

        dZ = dZ_new;
    }

    // Longitude, latitude, altitude
    lon = atan2 ( Y, X );
    lat = atan2 ( ZdZ, sqrt(rho2) );
    h   = Nh - N;
}