#ifndef PROYECTO_POSITION_H
#define PROYECTO_POSITION_H

#include "Matrix.h"
#include "SAT_const.h"
#include <cmath>

//------------------------------------------------------------------------------
// Matrix Position(double lon, double lat, double h)
//------------------------------------------------------------------------------
/**
 * Computes the position vector in the Earth-fixed frame.
 *
 * @param lon  Longitude (radians).
 * @param lat  Latitude (radians).
 * @param h    Height above the reference ellipsoid (meters).
 * @return     Position vector.
 */
//------------------------------------------------------------------------------
Matrix Position(double lon, double lat, double h);


#endif //PROYECTO_POSITION_H
