#ifndef PROYECTO_POLEMATRIX_H
#define PROYECTO_POLEMATRIX_H

#include "../include/Matrix.h"
#include "../include/R_x.h"
#include "../include/R_y.h"
#include "../include/R_z.h"

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
Matrix PoleMatrix (double xp, double yp);

#endif //PROYECTO_POLEMATRIX_H
