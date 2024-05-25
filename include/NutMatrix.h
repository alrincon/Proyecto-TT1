#ifndef PROYECTO_NUTMATRIX_H
#define PROYECTO_NUTMATRIX_H

#include "../include/Matrix.h"
#include "../include/R_x.h"
#include "../include/R_y.h"
#include "../include/R_z.h"
#include "../include/MeanObliquity.h"
#include "../include/NutAngles.h"

//------------------------------------------------------------------------------
// Matrix NutMatrix(double Mjd_TT)
//------------------------------------------------------------------------------
/**
 * Computes the transformation matrix from the mean to the true equator and equinox.
 *
 * @param Mjd_TT  Modified Julian Date (MJD) in Terrestrial Time (TT).
 * @return        Transformation matrix.
 */
//------------------------------------------------------------------------------
Matrix NutMatrix (double Mjd_TT);

#endif
