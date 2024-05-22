#include "../include/dotProduct.h"

double dotProduct(Matrix* v1, Matrix* v2){
    double dot;
    if(v1->getFilas() == 1 && v2->getFilas() == 1){
        dot = (*v1)(1,1)*(*v2)(1,1)+(*v1)(1,2)*(*v2)(1,2)+(*v1)(1,3)*(*v2)(1,3);
    }

    if(v1->getColumnas() == 1 && v2->getColumnas() == 1){
        dot = (*v1)(1,1)*(*v2)(1,1)+(*v1)(2,1)*(*v2)(2,1)+(*v1)(3,1)*(*v2)(3,1);
    }

    if(v1->getFilas() == 1 && v2->getColumnas() == 1){
        dot = (*v1)(1,1)*(*v2)(1,1)+(*v1)(1,2)*(*v2)(2,1)+(*v1)(1,3)*(*v2)(3,1);
    }

    if(v1->getColumnas() == 1 && v2->getFilas() == 1){
        dot = (*v1)(1,1)*(*v2)(1,1)+(*v1)(2,1)*(*v2)(1,2)+(*v1)(3,1)*(*v2)(1,3);
    }

    //Matrix v(1,1);
    //v = (*v1)*((*v2).transpose());
    //return(v(1,1));
    return dot;
}