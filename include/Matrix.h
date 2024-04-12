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
    double& operator()(const int i, const int j) const;

    void print();

private:
    void initMatrix();

private:
    int fil;
    int col;
    double **matrix;
};

#endif

