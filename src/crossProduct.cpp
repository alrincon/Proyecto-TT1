#include "../include/crossProduct.h"

//entrada (1,3) (1,3) salida (1,3)
Matrix crossProduct(Matrix* v1, Matrix* v2){

    if(v1->getFilas() != 1 || v2->getFilas() != 1 || v1->getColumnas() != 3 || v2->getColumnas() != 3){
        cout << "Error, dimensiones paar crossProduct incorrectas" << endl;
    }

    Matrix result(1, 3);
    result(1,1) = (*v1)(1,2)*(*v2)(1,3) - (*v1)(1,3)*(*v2)(1,2);
    result(1,2) = (*v1)(1,3)*(*v2)(1,1) - (*v1)(1,1)*(*v2)(1,3);
    result(1,3) = (*v1)(1,1)*(*v2)(1,2) - (*v1)(1,2)*(*v2)(1,1);

    return result;
}