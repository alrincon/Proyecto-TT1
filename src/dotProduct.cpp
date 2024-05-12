#include "../include/dotProduct.h"

double dotProduct(Matrix* v1, Matrix* v2){
    Matrix v = (*v1)*(*v2).transpose();
    return(v(1,1));
}