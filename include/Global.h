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
int infFile = 100;

Matrix eopdata(infFile, 13);
Matrix Cnm(181, 181);
Matrix Snm(181, 181);
Matrix PC(1,1);
Matrix obs(nobs, 4);


#endif
