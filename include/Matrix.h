#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <stdexcept>

class Matrix {
private:
    std::vector<std::vector<double>> datos;
    int filas;
    int columnas;

public:
    // Constructor
    Matrix(int filas, int columnas);

    // Métodos para obtener el número de filas y columnas
    int numColumnas() const;
    int numFilas() const;

    // Métodos para establecer y obtener valores
    void set(int row, int column, double value);
    double get(int row, int column) const;

    // Operadores aritméticos
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(double scalar) const;
    Matrix transpose() const;
    Matrix operator*(const Matrix& other) const;

    // Método para resolver sistemas de ecuaciones lineales
    std::vector<double> solve(std::vector<double> b) const;
};

#endif // MATRIX_H
