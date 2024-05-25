#include "../include/angl.h"
#include "../include/dotProduct.h"

//------------------------------------------------------------------------------
// double sign(double x)
//------------------------------------------------------------------------------
/*
 * Returns the sign of a number.
 * This function returns -1 if the input is negative, otherwise it returns 1.
 *
 * @param <x> the input value
 * @return -1 if the input is negative, otherwise 1
 * @exception none
 * @note useful for ensuring the correct sign in calculations
 */
//------------------------------------------------------------------------------
double sign (double x){
    if(x < 0){
        return -x;
    }

    return x;
}

//------------------------------------------------------------------------------
// double angl(Matrix* vec1, Matrix* vec2)
//------------------------------------------------------------------------------
/**
 * Calculates the angle between two vectors.
 * This function computes the angle between two vectors in radians using the dot product.
 *
 * @param <vec1> pointer to the first vector (Matrix)
 * @param <vec2> pointer to the second vector (Matrix)
 * @return the angle between the two vectors in radians, or a large value if undefined
 * @exception none
 * @note returns a large value if either vector is zero-length
 */
//------------------------------------------------------------------------------
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