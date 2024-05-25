#include "../include/unit.h"

#include <stdio.h>
#include <vector>
#include <iostream>

/**
 * @brief Normalizes a vector to have unit length, unless it is near zero in which case a zero vector is returned.
 *
 * This function calculates the unit vector of the provided vector by dividing it by its norm (magnitude).
 * If the magnitude of the vector is very small (less than a predefined small threshold), the function returns
 * a zero vector to avoid division by zero and to handle near-zero vectors appropriately.
 *
 * @param vec Pointer to the Matrix object representing the vector to be normalized.
 * @return Matrix A Matrix object that either contains the unit vector or is a zero vector.
 *
 * @note The function assumes that the input Matrix is a vector (either a row or column vector).
 *       The function uses a predefined threshold ('small') of 0.000001 to determine if the vector is
 *       effectively zero. Adjust this threshold based on the precision requirements of your application.
 */
Matrix unit (Matrix* vec){
    double small = 0.000001;
    double magv = (*vec).norm();

    if ( magv > small ) {
        return (*vec)*(1/magv);
    }else {
        return (*vec)*0;
    }
}
