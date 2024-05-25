#include "../include/crossProduct.h"

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
Matrix crossProduct(Matrix* v1, Matrix* v2){
    Matrix result(v1->getFilas(), v1->getColumnas());

    if(v1->getFilas() == 1 && v1->getColumnas() == 3 && v2->getFilas() == 1 && v2->getColumnas() == 3 ){
        result(1,1) = (*v1)(1,2)*(*v2)(1,3) - (*v1)(1,3)*(*v2)(1,2);
        result(1,2) = (*v1)(1,3)*(*v2)(1,1) - (*v1)(1,1)*(*v2)(1,3);
        result(1,3) = (*v1)(1,1)*(*v2)(1,2) - (*v1)(1,2)*(*v2)(1,1);
    }else{
        if(v1->getFilas() == 3 && v1->getColumnas() == 1 && v2->getFilas() == 3 && v2->getColumnas() == 1 ){
            result(1,1) = (*v1)(2,1)*(*v2)(3,1) - (*v1)(3,1)*(*v2)(2,1);
            result(2,1) = (*v1)(3,1)*(*v2)(1,1) - (*v1)(1,1)*(*v2)(3,1);
            result(3,1) = (*v1)(1,1)*(*v2)(2,1) - (*v1)(2,1)*(*v2)(1,1);
        }else{
            cout << "Error, dimensiones paar crossProduct incorrectas" << endl;
        }
    }

    return result;
}