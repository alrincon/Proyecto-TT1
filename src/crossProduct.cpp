#include "../include/crossProduct.h"


Matrix crossProduct(Matrix* v1, Matrix* v2){
    Matrix result(1, 3);
    result(1,1) = (*v1)(1,2)*(*v2)(1,3) - (*v1)(1,3)*(*v2)(1,2);
    result(1,2) = (*v1)(1,3)*(*v2)(1,1) - (*v1)(1,1)*(*v2)(1,3);
    result(1,3) = (*v1)(1,1)*(*v2)(1,2) - (*v1)(1,2)*(*v2)(1,1);

    return result;
}