#ifndef SAT_CONST_H
#define SAT_CONST_H

namespace SAT_Const {
    // Mathematical constants
    extern const double pi2;
    extern const double Rad;
    extern const double Deg;
    extern const double Arcs;

    // General
    extern const double MJD_J2000;
    extern const double T_B1950;
    extern const double c_light;
    extern const double AU;

    // Physical parameters of the Earth, Sun and Moon
    extern const double R_Earth;
    extern const double f_Earth;
    extern const double R_Sun;
    extern const double R_Moon;

    // Earth rotation
    extern const double omega_Earth;

    // Gravitational coefficients
    extern const double GM_Earth;
    extern const double GM_Sun;
    extern const double GM_Moon;
    extern const double GM_Mercury;
    extern const double GM_Venus;
    extern const double GM_Mars;
    extern const double GM_Jupiter;
    extern const double GM_Saturn;
    extern const double GM_Uranus;
    extern const double GM_Neptune;
    extern const double GM_Pluto;

    // Solar radiation pressure at 1 AU
    extern const double P_Sol;
}

#endif
