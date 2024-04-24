#include "../include/sign_.h"
#include <cmath>

double sign(double a, double b){
    if(b>=0.0){
        return abs(a);
    } else{
        return -abs(a);
    }
}
