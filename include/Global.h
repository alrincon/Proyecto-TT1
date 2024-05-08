#ifndef PROYECTO_GLOBAL_H
#define PROYECTO_GLOBAL_H

#include "../include/Matrix.h"

Matrix eopdata(13, 6);
Matrix GEOS3data(13, 6);

Matrix Cnm(181, 181);
Matrix Snm(181, 181);

struct AUXPARAM{
    char *Mjd_UTD;
    int n;
    int m;
    int sun;
    int moon;
    int planets;
};


#endif
