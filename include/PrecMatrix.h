#ifndef PROYECTO_PRECMATRIX_H
#define PROYECTO_PRECMATRIX_H

#include "../include/Matrix.h"
#include "../include/SAT_const.h"
#include "../include/R_x.h"
#include "../include/R_y.h"
#include "../include/R_z.h"

//------------------------------------------------------------------------------
// Matrix PrecMatrix(double Mjd_1, double Mjd_2)
//------------------------------------------------------------------------------
/**
 * Computes the precession transformation matrix.
 *
 * @param Mjd_1  Modified Julian Date (MJD) of initial epoch.
 * @param Mjd_2  Modified Julian Date (MJD) of final epoch.
 * @return       Precession transformation matrix.
 */
//------------------------------------------------------------------------------
Matrix  PrecMatrix (double Mjd_1, double Mjd_2);


#endif //PROYECTO_PRECMATRIX_H
