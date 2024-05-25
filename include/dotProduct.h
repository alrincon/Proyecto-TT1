#ifndef PROYECTO_DOTPRODUCT_H
#define PROYECTO_DOTPRODUCT_H

#include "Matrix.h"

//------------------------------------------------------------------------------
// dotProduct(Matrix* v1, Matrix* v2)
//------------------------------------------------------------------------------
/**
 * Computes the dot product of two vectors.
 *
 * This function calculates the dot product of two vectors. It accepts two matrices
 * representing vectors and returns their dot product.
 *
 * @param v1 Pointer to the first vector matrix.
 * @param v2 Pointer to the second vector matrix.
 * @return Dot product of v1 and v2.
 *
 * @note Both input matrices must represent vectors of the same dimension.
 */
//------------------------------------------------------------------------------
double dotProduct(Matrix* v1, Matrix* v2);

#endif
