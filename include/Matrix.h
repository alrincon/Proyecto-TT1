#ifndef _MATRIX_
#define _MATRIX_

class Matrix
{
public:
    Matrix(int fil, int col);
    Matrix(int fil, int col, double v[], int n);
    Matrix(const Matrix& m);
    ~Matrix();

    Matrix& operator=(const Matrix& matrix2);
    Matrix  operator+(const Matrix& matrix2);
    Matrix  operator-(const Matrix& matrix2);
    Matrix  operator*(const Matrix& matrix2);
    Matrix  operator*(const double scalar);
    double& operator()(const int i, const int j) const;

    void print();

    int getFilas();
    int getColumnas();
    double norm();

private:
    void initMatrix();

protected:
    int fil;
    int col;
    double **matrix;
};

#endif

