#include "../include/PoleMatrix.h"

//------------------------------------------------------------------------------
// Matrix PoleMatrix(double xp, double yp)
//------------------------------------------------------------------------------
/**
 * Computes the transformation matrix for polar motion.
 *
 * @param xp  Polar motion coordinate in x-direction (arcseconds).
 * @param yp  Polar motion coordinate in y-direction (arcseconds).
 * @return    Transformation matrix.
 */
//------------------------------------------------------------------------------
Matrix PoleMatrix (double xp, double yp){
    return R_y(-xp) * R_x(-yp);
}