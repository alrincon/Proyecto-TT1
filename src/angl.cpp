#include "../include/angl.h"
#include "../include/dotProduct.h"

double sign (double x){
    if(x < 0){
        return -x;
    }

    return x;
}

double angl (Matrix* vec1, Matrix* vec2){
    double theta;

    double small     = 0.00000001;
    double undefined = 999999.1;

    double magv1 = (*vec1).norm();
    double magv2 = (*vec2).norm();

    double temp;

    if (magv1*magv2 > pow(small,2)) {
        temp = dotProduct(vec1, vec2) / (magv1 * magv2);


        if (abs(temp) > 1.0) {
            temp = sign(temp) * 1.0;
        }
        theta = acos(temp);
    }else {
        theta = undefined;
    }

    return theta;
}