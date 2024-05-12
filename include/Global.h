#ifndef PROYECTO_GLOBAL_H
#define PROYECTO_GLOBAL_H

#include "../include/Matrix.h"

int nobs = 46;
int infFile = 100;

Matrix eopdata(infFile, 13);

Matrix obs(nobs, 4);
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
