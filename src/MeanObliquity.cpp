#include "../include/MeanObliquity.h"

//------------------------------------------------------------------------------
// double MeanObliquity (double Mjd_TT)
//------------------------------------------------------------------------------
/**
 * Computes the mean obliquity of the ecliptic for a given Terrestrial Time (TT) Modified Julian Date.
 *
 * @param Mjd_TT Terrestrial Time (TT) Modified Julian Date.
 * @return       The mean obliquity of the ecliptic in radians.
 */
//------------------------------------------------------------------------------
double MeanObliquity (double Mjd_TT){
    double T = (Mjd_TT - MJD_J2000)/36525.0;

    double MOblq = Rad *( 84381.448/3600.0-(46.8150+(0.00059-0.001813*T)*T)*T/3600.0 );
    return MOblq;
}