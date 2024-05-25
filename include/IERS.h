#ifndef PROYECTO_IERS_H
#define PROYECTO_IERS_H


#include <iostream>
#include <vector>
#include "Matrix.h"
#include "SAT_const.h"

int findMatchRow(Matrix a, int b);

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
void IERS(Matrix* eop, double Mjd_UTC, char interp, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC);

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
void IERS(Matrix* eop, double Mjd_UTC, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC);


#endif //PROYECTO_IERS_H
