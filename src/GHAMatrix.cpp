#include "../include/GHAMatrix.h"


//------------------------------------------------------------------------------
// GHAMatrix (double Mjd_UT1)
//------------------------------------------------------------------------------
/**
 * Computes the Greenwich Hour Angle (GHA) matrix for a given Modified Julian
 * Date in UT1 time scale.
 *
 * @param Mjd_UT1 Modified Julian Date in UT1 time scale.
 * @return GHA matrix.
 */
//------------------------------------------------------------------------------
Matrix GHAMatrix (double Mjd_UT1){
    return R_z( gast(Mjd_UT1) );
}