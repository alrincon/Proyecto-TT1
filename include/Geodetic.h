#ifndef PROYECTO_GEODETIC_H
#define PROYECTO_GEODETIC_H

#include "Matrix.h"
#include <stdexcept>

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
void Geodetic(double& lon, double& lat, double& h, Matrix* r);

#endif //PROYECTO_GEODETIC_H
