#include "../include/IERS.h"
#include <stdio.h>
#include <vector>
#include <iostream>

//------------------------------------------------------------------------------
// void IERS(Matrix* eop, double Mjd_UTC, char interp,
//           double& x_pole, double& y_pole, double& UT1_UTC,
//           double& LOD, double& dpsi, double& deps,
//           double& dx_pole, double& dy_pole, double& TAI_UTC)
//------------------------------------------------------------------------------
/**
 * Retrieves Earth rotation parameters from the provided Earth Orientation
 * Parameters (EOP) data set for a given Modified Julian Date (Mjd_UTC) using
 * either linear interpolation or nearest neighbor interpolation.
 *
 * @param eop      Pointer to the Earth Orientation Parameters (EOP) data set.
 * @param Mjd_UTC  Modified Julian Date (MJD) in Coordinated Universal Time (UTC).
 * @param interp   Interpolation method ('l' for linear, 'n' for nearest neighbor).
 * @param x_pole   Pole coordinate x (radians).
 * @param y_pole   Pole coordinate y (radians).
 * @param UT1_UTC  UT1 - UTC time difference (seconds).
 * @param LOD      Length of day (seconds).
 * @param dpsi     Nutation correction in longitude (radians).
 * @param deps     Nutation correction in obliquity (radians).
 * @param dx_pole  Pole coordinate x rate (radians/day).
 * @param dy_pole  Pole coordinate y rate (radians/day).
 * @param TAI_UTC  TAI - UTC time difference (seconds).
 */
//------------------------------------------------------------------------------
void IERS(Matrix* eop, double Mjd_UTC, char interp, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC) {

    if (interp == 'l') {
        // linear interpolation
        double mjd = (floor(Mjd_UTC));


        int i = findMatchRow((*eop).extractRow(4), mjd);

        Matrix preeop = (*eop).extractCol(i);
        Matrix nexteop = (*eop).extractCol(i+1);

        double mfme = 1440.0 * (Mjd_UTC - floor(Mjd_UTC));
        double fixf = mfme / 1440.0;

        // Setting of IERS Earth rotation parameters  (UT1 - UTC[s], TAI - UTC[s], x["], y ["])
        x_pole = preeop(5,1) + (nexteop(5,1) - preeop(5,1)) * fixf;
        y_pole = preeop(6,1) + (nexteop(6,1) - preeop(6,1)) * fixf;
        UT1_UTC = preeop(7,1) + (nexteop(7,1) - preeop(7,1)) * fixf;
        LOD = preeop(8,1) + (nexteop(8,1) - preeop(8,1)) * fixf;
        dpsi = preeop(9,1) + (nexteop(9,1) - preeop(9,1)) * fixf;
        deps = preeop(10,1) + (nexteop(10,1) - preeop(10,1)) * fixf;
        dx_pole = preeop(11,1) + (nexteop(11,1)- preeop(11,1)) * fixf;
        dy_pole = preeop(12,1) + (nexteop(12,1) - preeop(12,1)) * fixf;
        TAI_UTC = preeop(13,1);

        x_pole = x_pole /Arcs;  // Pole coordinate[rad]
        y_pole = y_pole /Arcs;  // Pole coordinate[rad]
        dpsi = dpsi /Arcs;
        deps = deps /Arcs;
        dx_pole = dx_pole /Arcs; // Pole coordinate[rad]
        dy_pole = dy_pole /Arcs; // Pole coordinate[rad]
    }else if (interp =='n'){
        double mjd = (floor(Mjd_UTC));

        int i = findMatchRow((*eop).extractRow(4), mjd);

        Matrix teop = (*eop).extractCol(i);

        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])

        x_pole  = teop(5,1)/Arcs;  // Pole coordinate [rad]
        y_pole  = teop(6,1)/Arcs;  // Pole coordinate [rad]
        UT1_UTC = teop(7,1);             // UT1-UTC time difference [s]
        LOD     = teop(8,1);             // Length of day [s]
        dpsi    = teop(9,1)/Arcs;
        deps    = teop(10,1)/Arcs;
        dx_pole = teop(11,1)/Arcs; // Pole coordinate [rad]
        dy_pole = teop(12,1)/Arcs; // Pole coordinate [rad]
        TAI_UTC = teop(13,1);            // TAI-UTC time difference [s]
    }
}

//------------------------------------------------------------------------------
// void IERS(Matrix* eop, double Mjd_UTC, double& x_pole, double& y_pole,
//           double& UT1_UTC, double& LOD, double& dpsi, double& deps,
//           double& dx_pole, double& dy_pole, double& TAI_UTC)
//------------------------------------------------------------------------------
/**
 * Retrieves Earth rotation parameters from the provided Earth Orientation
 * Parameters (EOP) data set for a given Modified Julian Date (Mjd_UTC) using
 * nearest neighbor interpolation by default.
 *
 * @param eop      Pointer to the Earth Orientation Parameters (EOP) data set.
 * @param Mjd_UTC  Modified Julian Date (MJD) in Coordinated Universal Time (UTC).
 * @param x_pole   Pole coordinate x (radians).
 * @param y_pole   Pole coordinate y (radians).
 * @param UT1_UTC  UT1 - UTC time difference (seconds).
 * @param LOD      Length of day (seconds).
 * @param dpsi     Nutation correction in longitude (radians).
 * @param deps     Nutation correction in obliquity (radians).
 * @param dx_pole  Pole coordinate x rate (radians/day).
 * @param dy_pole  Pole coordinate y rate (radians/day).
 * @param TAI_UTC  TAI - UTC time difference (seconds).
 */
//------------------------------------------------------------------------------
void IERS(Matrix* eop, double Mjd_UTC, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC){
    return IERS(eop, Mjd_UTC, 'n', x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
}

/**
 * Finds the index of the row in a matrix `a` that matches the integer `b`.
 *
 * @param a The matrix to search.
 * @param b The integer value to match.
 * @return The index of the matching row, or 0 if no match is found.
 */
//------------------------------------------------------------------------------
int findMatchRow(Matrix a, int b){
    for(int i = 1; i <= a.getColumnas(); i++){
        if(abs(a(1,i) - b) < 10e-6){
            return i;
        }
    }

    return 0;
}


