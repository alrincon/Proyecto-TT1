#ifndef _MATRIX_
#define _MATRIX_

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

class Matrix
{
public:
    //------------------------------------------------------------------------------
// Matrix::Matrix(int fil, int col)
//------------------------------------------------------------------------------
/**
 * Constructs a matrix with the specified number of rows and columns, initialized to zeros.
 *
 * @param fil Number of rows.
 * @param col Number of columns.
 */
//------------------------------------------------------------------------------
    Matrix(int fil, int col);
    //------------------------------------------------------------------------------
// Matrix::Matrix(int fil, int col, double v[], int n)
//------------------------------------------------------------------------------
/**
 * Constructs a matrix with the specified number of rows and columns, initialized with provided values. If the number of provided values is less than the total size of the matrix, the remaining elements are set to zero.
 *
 * @param fil Number of rows.
 * @param col Number of columns.
 * @param v   Array containing the values to initialize the matrix.
 * @param n   Number of values in the array.
 */
//------------------------------------------------------------------------------
    Matrix(int fil, int col, double v[], int n);
    //------------------------------------------------------------------------------
// Matrix::Matrix(int fil, int col, bool identity)
//------------------------------------------------------------------------------
/**
 * Constructs a matrix with the specified number of rows and columns, initialized as an identity matrix if 'identity' is true, or as a zero matrix otherwise.
 *
 * @param fil      Number of rows.
 * @param col      Number of columns.
 * @param identity If true, initializes the matrix as an identity matrix; otherwise, initializes it as a zero matrix.
 */
//------------------------------------------------------------------------------
    Matrix(int fil, int col, bool identity);
    //------------------------------------------------------------------------------
// Matrix::Matrix(const Matrix& m)
//------------------------------------------------------------------------------
/**
 * Constructs a matrix as a copy of another matrix.
 *
 * @param m The matrix to be copied.
 */
//------------------------------------------------------------------------------
    Matrix(const Matrix& m);
    //------------------------------------------------------------------------------
// Matrix::~Matrix()
//------------------------------------------------------------------------------
/**
 * Destructor for the Matrix class.
 */
//------------------------------------------------------------------------------
    ~Matrix();

    //------------------------------------------------------------------------------
// Matrix& Matrix::operator=(const Matrix& matrix2)
//------------------------------------------------------------------------------
/**
 * Assigns the values of another matrix to this matrix.
 *
 * @param matrix2 The matrix whose values are to be assigned.
 * @return        Reference to this matrix after assignment.
 */
//------------------------------------------------------------------------------
    Matrix& operator=(const Matrix& matrix2);
    //------------------------------------------------------------------------------
// Matrix Matrix::operator+(const Matrix& matrix2)
//------------------------------------------------------------------------------
/**
 * Adds another matrix to this matrix.
 *
 * @param matrix2 The matrix to be added.
 * @return        The result of the addition.
 */
//------------------------------------------------------------------------------
    Matrix  operator+(const Matrix& matrix2);
    //------------------------------------------------------------------------------
// Matrix Matrix::operator-(const Matrix& matrix2)
//------------------------------------------------------------------------------
/**
 * Subtracts another matrix from this matrix.
 *
 * @param matrix2 The matrix to be subtracted.
 * @return        The result of the subtraction.
 */
//------------------------------------------------------------------------------
    Matrix  operator-(const Matrix& matrix2);
    //------------------------------------------------------------------------------
// Matrix Matrix::operator*(const Matrix& matrix2)
//------------------------------------------------------------------------------
/**
 * Multiplies this matrix by another matrix.
 *
 * @param matrix2 The matrix to be multiplied.
 * @return        The result of the multiplication.
 */
//------------------------------------------------------------------------------
    Matrix  operator*(const Matrix& matrix2);
    //------------------------------------------------------------------------------
// Matrix Matrix::operator*(const double scalar)
//------------------------------------------------------------------------------
/**
 * Multiplies each element of the matrix by a scalar value.
 *
 * @param scalar The scalar value to multiply by.
 * @return       The result of the multiplication.
 */
//------------------------------------------------------------------------------
    Matrix  operator*(const double scalar);
    //------------------------------------------------------------------------------
// double& Matrix::operator()(const int i, const int j) const
//------------------------------------------------------------------------------
/**
 * Accesses the element at the specified row and column of the matrix.
 *
 * @param i Row index (1-based).
 * @param j Column index (1-based).
 * @return  Reference to the element at the specified row and column.
 */
//------------------------------------------------------------------------------
    double& operator()(const int i, const int j) const;
//------------------------------------------------------------------------------
// void Matrix::setTam(int nfil, int ncol)
//------------------------------------------------------------------------------
/**
 * Sets the number of rows and columns in the matrix and initializes it with zeros.
 *
 * @param nfil Number of rows.
 * @param ncol Number of columns.
 */
//------------------------------------------------------------------------------
    void setTam(int nfil, int ncol);
//------------------------------------------------------------------------------
// void Matrix::redefine(Matrix* mat)
//------------------------------------------------------------------------------
/**
 * Redefines the matrix using the values of another matrix.
 *
 * @param mat The matrix whose values are to be copied.
 */
//------------------------------------------------------------------------------
    void redefine(Matrix* mat);

    //------------------------------------------------------------------------------
// void Matrix::print()
//------------------------------------------------------------------------------
/**
 * Prints the matrix to the standard output.
 */
//------------------------------------------------------------------------------
    void print();
//------------------------------------------------------------------------------
// int Matrix::getFilas()
//------------------------------------------------------------------------------
/**
 * Gets the number of rows in the matrix.
 *
 * @return Number of rows.
 */
//------------------------------------------------------------------------------
    int getFilas();
    //------------------------------------------------------------------------------
// int Matrix::getColumnas()
//------------------------------------------------------------------------------
/**
 * Gets the number of columns in the matrix.
 *
 * @return Number of columns.
 */
//------------------------------------------------------------------------------
    int getColumnas();
    //------------------------------------------------------------------------------
// double Matrix::norm()
//------------------------------------------------------------------------------
/**
 * Computes the norm of the matrix.
 *
 * @return The norm of the matrix.
 */
//------------------------------------------------------------------------------
    double norm();
    //------------------------------------------------------------------------------
// void Matrix::assign(int i, int j, double v)
//------------------------------------------------------------------------------
/**
 * Assigns a value to the element at the specified row and column of the matrix.
 *
 * @param i Row index (1-based).
 * @param j Column index (1-based).
 * @param v Value to be assigned.
 */
//------------------------------------------------------------------------------
    void assign(int i, int j, double v);
//------------------------------------------------------------------------------
// Matrix Matrix::extractCol(int i)
//------------------------------------------------------------------------------
/**
 * Extracts a column from the matrix as a new matrix.
 *
 * @param i Index of the column to be extracted (1-based).
 * @return  A matrix containing the extracted column.
 */
//------------------------------------------------------------------------------
    Matrix extractCol(int i);
    //------------------------------------------------------------------------------
// Matrix Matrix::extractRow(int i)
//------------------------------------------------------------------------------
/**
 * Extracts a row from the matrix as a new matrix.
 *
 * @param i Index of the row to be extracted (1-based).
 * @return  A matrix containing the extracted row.
 */
//------------------------------------------------------------------------------
    Matrix extractRow(int i);
    //------------------------------------------------------------------------------
// Matrix Matrix::transpose()
//------------------------------------------------------------------------------
/**
 * Computes the transpose of the matrix.
 *
 * @return The transpose of the matrix.
 */
//------------------------------------------------------------------------------
    Matrix transpose();
    static double dotProduct(Matrix* v1, Matrix* v2);
    //------------------------------------------------------------------------------
// Matrix Matrix::inverse()
//------------------------------------------------------------------------------
/**
 * Computes the inverse of the matrix.
 *
 * @return The inverse of the matrix.
 */
//------------------------------------------------------------------------------
    Matrix inverse();
    //------------------------------------------------------------------------------
// Matrix Matrix::clone()
//------------------------------------------------------------------------------
/**
 * Creates a deep copy of the matrix.
 *
 * @return A copy of the matrix.
 */
//------------------------------------------------------------------------------
    Matrix clone();

private:
    void initMatrix();

protected:
    int fil;
    int col;
    double **matrix;
};

#endif

