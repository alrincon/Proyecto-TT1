using namespace std;
#include "Matrix.h"

#include <vector>
#include <stdexcept>

Matrix::Matrix(int filas, int columnas) : filas(filas), columnas(columnas) {
    datos.resize(filas, vector<double>(columnas, 0.0));
}

// Método para get el número de filas
int Matrix::numFilas() const {
    return filas;
}

// Método para get el número de columnas
int Matrix::numColumnas() const {
    return columnas;
}

// Método para set un valor a un elemento de la Matrix
void Matrix::set(int fila, int columna, double valor) {
    datos[fila][columna] = valor;
}

// Método para get un elemento de la Matrix
double Matrix::get(int fila, int columna) const {
    return datos[fila][columna];
}

// Operador de suma
Matrix Matrix::operator+(const Matrix& otra) const {
    if (filas != otra.filas || columnas != otra.columnas) {
        throw invalid_argument("Las dimensiones de las matrices no coinciden.");
    }

    Matrix resultado(filas, columnas);
    for (int i = 0; i < filas; ++i) {
        for (int j = 0; j < columnas; ++j) {
            resultado.datos[i][j] = datos[i][j] + otra.datos[i][j];
        }
    }
    return resultado;
}

// Operador de resta
Matrix Matrix::operator-(const Matrix& otra) const {
    if (filas != otra.filas || columnas != otra.columnas) {
        throw invalid_argument("Las dimensiones de las matrices no coinciden.");
    }

    Matrix resultado(filas, columnas);
    for (int i = 0; i < filas; ++i) {
        for (int j = 0; j < columnas; ++j) {
            resultado.datos[i][j] = datos[i][j] - otra.datos[i][j];
        }
    }
    return resultado;
}

// Operador de multiplicación por un escalar
Matrix Matrix::operator*(double escalar) const {
    Matrix resultado(filas, columnas);
    for (int i = 0; i < filas; ++i) {
        for (int j = 0; j < columnas; ++j) {
            resultado.datos[i][j] = datos[i][j] * escalar;
        }
    }
    return resultado;
}

// Operador de transposición
Matrix Matrix::transpose() const {
    Matrix resultado(columnas, filas);
    for (int i = 0; i < filas; ++i) {
        for (int j = 0; j < columnas; ++j) {
            resultado.datos[j][i] = datos[i][j];
        }
    }
    return resultado;
}

// Operador de multiplicación de matrices
Matrix Matrix::operator*(const Matrix& otra) const {
    if (columnas != otra.filas) {
        throw invalid_argument("El número de columnas de la primera Matrix debe ser igual al número de filas de la segunda Matrix.");
    }

    Matrix resultado(filas, otra.columnas);
    for (int i = 0; i < filas; ++i) {
        for (int j = 0; j < otra.columnas; ++j) {
            double sum = 0.0;
            for (int k = 0; k < columnas; ++k) {
                sum += datos[i][k] * otra.datos[k][j];
            }
            resultado.datos[i][j] = sum;
        }
    }
    return resultado;
}

// Resolver sistema de ecuaciones lineales
vector<double> Matrix::solve(vector<double> b) const {
    if (filas != columnas || filas != b.size()) {
        throw invalid_argument("La Matrix debe ser cuadrada y tener el mismo número de filas que el vector b.");
    }

    // Implementación básica: usar la eliminación de Gauss sin pivoteo
    Matrix A_extendida(filas, columnas + 1);
    for (int i = 0; i < filas; ++i) {
        for (int j = 0; j < columnas; ++j) {
            A_extendida.set(i, j, datos[i][j]);
        }
        A_extendida.set(i, columnas, b[i]);
    }

    for (int i = 0; i < filas; ++i) {
        // Hacer que el pivote sea 1
        double pivot = A_extendida.get(i, i);
        if (pivot == 0) {
            throw runtime_error("El pivote es cero. Se necesita pivoteo para resolver el sistema.");
        }
        for (int j = i; j <= columnas; ++j) {
            A_extendida.set(i, j, A_extendida.get(i, j) / pivot);
        }

        // Eliminación hacia adelante
        for (int k = i + 1; k < filas; ++k) {
            double factor = A_extendida.get(k, i);
            for (int j = i; j <= columnas; ++j) {
                A_extendida.set(k, j, A_extendida.get(k, j) - factor * A_extendida.get(i, j));
            }
        }
    }

    // Sustitución hacia atrás
    vector<double> x(filas);
    for (int i = filas - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < columnas; ++j) {
            sum += A_extendida.get(i, j) * x[j];
        }
        x[i] = A_extendida.get(i, columnas) - sum;
    }

    return x;
}
