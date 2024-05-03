#include <cmath>

namespace SAT_Const {
    // Mathematical constants
    const double pi2 = 2 * M_PI;
    const double Rad = M_PI / 180;
    const double Deg = 180 / M_PI;
    const double Arcs = 3600 * 180 / M_PI;

    // General
    const double MJD_J2000 = 51544.5;
    const double T_B1950 = -0.500002108;
    const double c_light = 299792458.000000000;
    const double AU = 149597870700.000000;

    // Physical parameters of the Earth, Sun and Moon
    const double R_Earth = 6378.1363e3;
    const double f_Earth = 1 / 298.257223563;
    const double R_Sun = 696000e3;
    const double R_Moon = 1738e3;

    // Earth rotation
    const double omega_Earth = 15.04106717866910 / 3600 * Rad;

    // Gravitational coefficients
    const double GM_Earth = 398600.435436e9;
    const double GM_Sun = 132712440041.939400e9;
    const double GM_Moon = GM_Earth / 81.30056907419062;
    const double GM_Mercury = 22031.780000e9;
    const double GM_Venus = 324858.592000e9;
    const double GM_Mars = 42828.375214e9;
    const double GM_Jupiter = 126712764.800000e9;
    const double GM_Saturn = 37940585.200000e9;
    const double GM_Uranus = 5794548.600000e9;
    const double GM_Neptune = 6836527.100580e9;
    const double GM_Pluto = 977.0000000000009e9;

    // Solar radiation pressure at 1 AU
    const double P_Sol = 1367 / c_light;
}
