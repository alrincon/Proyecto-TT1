#ifndef PROYECTO_LEGENDRE_H
#define PROYECTO_LEGENDRE_H
#include <iostream>
#include <vector>
#include "Matrix.h"

//------------------------------------------------------------------------------
// void Legendre(int n, int m, double fi, Matrix& pnm, Matrix& dpnm)
//------------------------------------------------------------------------------
/**
 * Computes the associated Legendre functions and their derivatives up to degree n and order m at a given latitude.
 *
 * @param n     Maximum degree of the Legendre functions.
 * @param m     Maximum order of the Legendre functions.
 * @param fi    Latitude in radians at which the Legendre functions are evaluated.
 * @param pnm   Matrix storing the Legendre functions up to degree n and order m.
 * @param dpnm  Matrix storing the derivatives of the Legendre functions up to degree n and order m.
 */
//------------------------------------------------------------------------------
void  Legendre(int n, int m, double fi, Matrix& pnm, Matrix& dpnm);


#endif //PROYECTO_LEGENDRE_H
