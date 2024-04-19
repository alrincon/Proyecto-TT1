#include <stdexcept>
#include "../include/Vector.h"

Vector::Vector(int size): Matrix(size, 1){
}

Vector::Vector(int size, double v[], int n): Matrix(size, 1, v, n){
}

Vector::Vector(const Vector& v): Matrix((Matrix) v){
}

Vector::Vector(Matrix matrix): Matrix(matrix) {
    if(matrix.getFilas() != 1){
        throw std::invalid_argument("La matriz debe ser un vector para convertirla");
    }
}

double Vector::norm(){
    double suma = 0.0;

    for(int i = 0; i < fil; i ++){
        suma += matrix[i][0]*matrix[i][0];
    }

    return sqrt(suma);
}

Vector& Vector::operator=(const Vector& vector2){
    Matrix::operator=(vector2);
    return *this;
}

Vector Vector::operator+(const Vector& vector2){
    return static_cast<Vector>(Matrix::operator+(vector2));
}

Vector Vector::operator-(const Vector& vector2){
    return static_cast<Vector>(Matrix::operator-(vector2));
}