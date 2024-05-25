#include "../include/LTC.h"

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
Matrix LTC(double lon, double lat){
    Matrix M = R_y(-1.0*lat)*R_z(lon);

    for (int j = 1; j <=3; j++){
        double Aux = M(1, j);

        M(1, j) = M(2, j);
        M(2, j) = M(3, j);
        M(3, j) = Aux;
    }

    return M;
}
