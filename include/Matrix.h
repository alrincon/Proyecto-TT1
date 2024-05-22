#ifndef _MATRIX_
#define _MATRIX_

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

class Matrix
{
public:
    Matrix(int fil, int col);
    Matrix(int fil, int col, double v[], int n);
    Matrix(int fil, int col, bool identity);
    Matrix(const Matrix& m);
    ~Matrix();

    Matrix& operator=(const Matrix& matrix2);
    Matrix  operator+(const Matrix& matrix2);
    Matrix  operator-(const Matrix& matrix2);
    Matrix  operator*(const Matrix& matrix2);
    Matrix  operator*(const double scalar);
    double& operator()(const int i, const int j) const;

    void setTam(int nfil, int ncol);

    void redefine(Matrix* mat);

    void print();

    int getFilas();
    int getColumnas();
    double norm();
    void assign(int i, int j, double v);

    Matrix extractCol(int i);
    Matrix extractRow(int i);
    Matrix transpose();
    static double dotProduct(Matrix* v1, Matrix* v2);
    Matrix inverse();
    Matrix clone();

private:
    void initMatrix();

protected:
    int fil;
    int col;
    double **matrix;
};

#endif

