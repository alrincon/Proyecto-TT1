#ifndef PROYECTO_GLOBAL_H
#define PROYECTO_GLOBAL_H

#include "../include/Matrix.h"
#include "types.h"
/*
class Global{
public:


    static void initialize();
};
*/

aux AuxParam;

int nobs = 46;
int infFile = 21413;

Matrix eopdata(13,infFile);
Matrix Cnm(181, 181);
Matrix Snm(181, 181);
Matrix PC(2285,1020);
Matrix obs(nobs, 4);


#endif
