#include "../include/unit.h"

#include <stdio.h>
#include <vector>
#include <iostream>

Matrix unit (Matrix* vec){
    double small = 0.000001;
    double magv = (*vec).norm();

    if ( magv > small ) {
        return (*vec)*(1/magv);
    }else {
        return (*vec)*0;
    }
}
