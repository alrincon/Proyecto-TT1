#include "../include/EqnEquinox.h"


//------------------------------------------------------------------------------
// EqnEquinox(double Mjd_TT)
//------------------------------------------------------------------------------
/**
 * Computes the equation of the equinoxes.
 *
 * This function computes the equation of the equinoxes, which is the
 * difference between apparent and mean sidereal times, using the nutation in
 * longitude and the mean obliquity of the ecliptic.
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time).
 * @return Equation of the equinoxes.
 */
//------------------------------------------------------------------------------
double EqnEquinox (double Mjd_TT){
    // Nutation in longitude and obliquity
    double dpsi;
    double deps;
    NutAngles (Mjd_TT, dpsi, deps);

    // Equation of the equinoxes
    double EqE = dpsi * cos ( MeanObliquity(Mjd_TT) );
    return EqE;
}