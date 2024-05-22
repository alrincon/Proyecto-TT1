#include "../include/realmod.h"

double realmod(double x, double y)
{
    if (y == 0 )
        return x;
    else if ( x == y)
        return 0.0;
    else
        return x - floor(x/y) *y;
}