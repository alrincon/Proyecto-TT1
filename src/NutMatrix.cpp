#include "../include/NutMatrix.h"

//------------------------------------------------------------------------------
// Matrix NutMatrix(double Mjd_TT)
//------------------------------------------------------------------------------
/**
 * Computes the transformation matrix from the mean to the true equator and equinox.
 *
 * @param Mjd_TT  Modified Julian Date (MJD) in Terrestrial Time (TT).
 * @return        Transformation matrix.
 */
//------------------------------------------------------------------------------
Matrix NutMatrix (double Mjd_TT){
    // Mean obliquity of the ecliptic
    double eps = MeanObliquity (Mjd_TT);

    // Nutation in longitude and obliquity
    double dpsi;
    double deps;
    NutAngles (Mjd_TT, dpsi, deps);

    // Transformation from mean to true equator and equinox
    return R_x(-eps-deps)*R_z(-dpsi)*R_x(+eps);
}