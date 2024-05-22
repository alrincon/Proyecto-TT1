#include "../include/Matrix.h"

Matrix::Matrix(int fil, int col) : fil(fil), col(col)
{
    initMatrix();
}

Matrix::Matrix(int fil, int col, double v[], int n): fil(fil), col(col)
{
    initMatrix();

    int k = 0;

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
            if (k < n)
                matrix[i][j] = v[k++];
            else
                matrix[i][j] = 0;
        }
}

Matrix::Matrix(int fil, int col, bool identity): fil(fil), col(col)
{
    if(identity){
        initMatrix();

        int k = 0;

        for (int i = 0; i < fil; i++){
            for (int j = 0; j < col; j++){
                if (i == j)
                    matrix[i][j] = 1;
                else
                    matrix[i][j] = 0;
            }
        }
    }else{
        initMatrix();

        int k = 0;

        for (int i = 0; i < fil; i++){
            for (int j = 0; j < col; j++){
                if (i == j)
                matrix[i][j] = 0;
                else
                matrix[i][j] = 0;
            }
        }
    }
}

Matrix::Matrix(const Matrix& m)
{
    *this = m;
}

Matrix::~Matrix()
{
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];

    delete[] matrix;
}

void Matrix::initMatrix()
{
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}

Matrix& Matrix::operator=(const Matrix& matrix2)
{
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
            this->matrix[i][j] = matrix2.matrix[i][j];
        }

    return *this;
}

Matrix Matrix::operator+(const Matrix& matrix2)
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];

    return result;
}

Matrix Matrix::operator-(const Matrix& matrix2)
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];

    return result;
}

Matrix Matrix::operator*(const Matrix& matrix2)
{
    Matrix result(fil, col);

    for (int i = 0; i < this->fil ; i++){
        for (int j = 0; j < matrix2.col; j++){
            result.matrix[i][j] = 0;
            for (int k = 0; k < this->col; k++){
                result.matrix[i][j] = result.matrix[i][j] + this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }

    return result;
}


double& Matrix::operator()(const int i, const int j) const
{
    return matrix[i-1][j-1];
}

void Matrix::print()
{
    for (int i = 0; i < fil; i++){
        for (int j = 0; j < col; j++){
            std::cout << std::fixed << std::setprecision(14) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int Matrix::getFilas(){
    return fil;
}

int Matrix::getColumnas(){
    return col;
}

double Matrix::norm(){
    double suma = 0;

    if(fil == 1){
        for(int i = 0; i < col; i++){
            suma += pow(matrix[0][i],2);
        }
        return sqrt(suma);
    }

    if(col == 1){
        for(int i = 0; i < fil; i++){
            suma += pow(matrix[i][0],2);
        }
        return sqrt(suma);
    }

    if(col != 1 && fil != 1){
        cout << "Filas: " << fil << endl;
        cout << "Columnas: " << col << endl;
        throw runtime_error("An error occurred, fallo en las dimensiones de la norma");
    }

    return -1;
}

void Matrix::assign(int i, int j, double v){
    matrix[i][j] = v;
}

Matrix Matrix::extractCol(int i){
    Matrix res(fil, 1);

    for(int k = 0; k < fil; k++){
        res(k+1, 1) = matrix[k][i-1];
    }

    return res;
}

Matrix Matrix::extractRow(int i){
    Matrix res(1, col);

    for(int k = 0; k < col; k++){
        res(1, k+1) = matrix[i-1][k];
    }

    return res;
}

Matrix Matrix::operator*(const double scalar) {
    Matrix result(fil, col);
    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            result(i + 1, j+ 1) = matrix[i][j] * scalar;
        }
    }
    return result;
}


Matrix Matrix::transpose() {
    Matrix transposed(col, fil);
    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            transposed.matrix[j][i] = matrix[i][j];
        }
    }
    return transposed;
}
#include <algorithm>

Matrix Matrix::inverse() {
    if (fil != col) {
        throw std::runtime_error("Inverse can only be calculated for square matrices.");
    }
    int n = fil;
    Matrix inv(n, n);  // Create an identity matrix for the inverse
    for (int i = 0; i < n; ++i) {
        inv.matrix[i][i] = 1.0;
    }

    // Using a copy of the original matrix to apply the operations
    Matrix temp = (*this).clone();

    // Gaussian elimination with partial pivoting
    for (int i = 0; i < n; ++i) {
        // Pivoting
        double maxEl = std::abs(temp.matrix[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (std::abs(temp.matrix[k][i]) > maxEl) {
                maxEl = std::abs(temp.matrix[k][i]);
                maxRow = k;
            }
        }

        // Swap the maximum row with the current row (column by column)
        if (maxRow != i) {
            std::swap(temp.matrix[maxRow], temp.matrix[i]);
            std::swap(inv.matrix[maxRow], inv.matrix[i]);
        }

        // Set to zero the lower part of the current column
        for (int k = i + 1; k < n; k++) {
            double c = -temp.matrix[k][i] / temp.matrix[i][i];
            for (int j = i; j < n; j++) {
                if (i == j) {
                    temp.matrix[k][j] = 0;
                } else {
                    temp.matrix[k][j] += c * temp.matrix[i][j];
                }
            }
            for (int j = 0; j < n; j++) {
                inv.matrix[k][j] += c * inv.matrix[i][j];
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    for (int i = n - 1; i >= 0; i--) {
        double c = temp.matrix[i][i];
        if (c == 0) {
            throw std::runtime_error("Matrix is singular and cannot be inverted.");
        }
        for (int j = 0; j < n; j++) {
            inv.matrix[i][j] = inv.matrix[i][j] / c;
        }
        temp.matrix[i][i] = 1;

        for (int k = i - 1; k >= 0; k--) {
            c = -temp.matrix[k][i];
            for (int j = 0; j < n; j++) {
                inv.matrix[k][j] += c * inv.matrix[i][j];
            }
            temp.matrix[k][i] = 0;
        }
    }

    return inv;
}

Matrix Matrix::clone(){
    Matrix result(fil, col);
    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            result(i + 1, j+ 1) = matrix[i][j];
        }
    }
    return result;
}
