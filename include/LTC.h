#ifndef PROYECTO_LTC_H
#define PROYECTO_LTC_H

#include "Matrix.h"
#include "R_y.h"
#include "R_z.h"

//------------------------------------------------------------------------------
// Matrix LTC(double lon, double lat)
//------------------------------------------------------------------------------
/**
 * Computes the Local Tangent Coordinate (LTC) transformation matrix for a given longitude and latitude.
 * LTC transformation converts geocentric coordinates to local tangent coordinates.
 *
 * @param lon Longitude in radians.
 * @param lat Latitude in radians.
 * @return    LTC transformation matrix.
 */
//------------------------------------------------------------------------------
Matrix LTC(double lon, double lat);


#endif
