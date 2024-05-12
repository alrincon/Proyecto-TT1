#ifndef PROYECTO_MEASUPDATE_H
#define PROYECTO_MEASUPDATE_H

#include "Matrix.h"

//z, s, x vectores
void MeasUpdate(Matrix* z, Matrix* g, Matrix* s, Matrix*G , int n, Matrix& K, Matrix& x, Matrix& P);


#endif
