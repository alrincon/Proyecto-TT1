//
// Created by alrincon on 19/04/2024.
//

#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <vector>
#include "Matrix.h"

class Vector: public Matrix{
public:

    Vector(int size);
    Vector(int size, double v[], int n);
    Vector(const Vector& v);
    Vector(Matrix matrix);

    double norm();


    Vector& operator=(const Vector& vector2);
    Vector operator+(const Vector& vector2);
    Vector operator-(const Vector& vector2);
};


#endif VECTOR_H
