#ifndef PROYECTO_TIMEUPDATE_H
#define PROYECTO_TIMEUPDATE_H

#include "Matrix.h"

void TimeUpdate(Matrix& P, Matrix* Phi, Matrix* Qdt);
void TimeUpdate(Matrix& P, Matrix* Phi);

#endif
