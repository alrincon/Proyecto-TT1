#ifndef PROYECTO_CROSSPRODUCT_H
#define PROYECTO_CROSSPRODUCT_H

#include "Matrix.h"

//------------------------------------------------------------------------------
// crossProduct(Matrix* v1, Matrix* v2)
//------------------------------------------------------------------------------
/**
 * Computes the cross product of two 3-dimensional vectors.
 *
 * This function calculates the cross product of two 3-dimensional vectors.
 * It accepts two matrices representing vectors and returns their cross product.
 *
 * @param v1 Pointer to the first vector matrix.
 * @param v2 Pointer to the second vector matrix.
 * @return Matrix representing the cross product of v1 and v2.
 *
 * @note Both input matrices must represent 3-dimensional vectors.
 */
//------------------------------------------------------------------------------
Matrix crossProduct(Matrix* v1, Matrix* v2);

#endif
