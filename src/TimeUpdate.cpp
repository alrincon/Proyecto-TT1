#include "../include/TimeUpdate.h"

void TimeUpdate(Matrix& P, Matrix* Phi, Matrix* Qdt){
    P = (*Phi)*P*(*Phi).transpose() + *Qdt;
}


void TimeUpdate(Matrix& P, Matrix* Phi){
    P = (*Phi)*P*(*Phi).transpose();
}