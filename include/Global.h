#ifndef PROYECTO_GLOBAL_H
#define PROYECTO_GLOBAL_H

#include "../include/Matrix.h"

typedef struct{
    double Mjd_UTC;
    double Mjd_TT;
    int n;
    int m;
    int sun;
    int moon;
    int planets;
} aux;

class Global{
public:
    static aux AuxParam;

    static Matrix *eopdata;
    static Matrix *Cnm;
    static Matrix *Snm;
    static Matrix *PC;
    static Matrix *obs;

    static void initialize();
};

#endif
