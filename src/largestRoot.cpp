#include "../include/largestRoot.h"

#include <iostream>
#include <complex>
#include <vector>
#include <limits>
#include <cmath>

// Define a tolerance for checking if a number is real
const double TOL = 1e-10;

//------------------------------------------------------------------------------
// std::complex<double> evaluatePolynomial(const std::vector<std::complex<double>>& coeffs, const std::complex<double>& z)
//------------------------------------------------------------------------------
/**
 * Evaluates a polynomial with given coefficients at a specified complex point using Horner's method.
 *
 * @param coeffs Coefficients of the polynomial in descending order of degree.
 * @param z      Complex point at which the polynomial is evaluated.
 * @return       The value of the polynomial at the specified complex point.
 */
//------------------------------------------------------------------------------
std::complex<double> evaluatePolynomial(const std::vector<std::complex<double>>& coeffs, const std::complex<double>& z) {
    std::complex<double> result = coeffs[0];
    std::complex<double> power = 1.0;
    for (size_t i = 1; i < coeffs.size(); ++i) {
        power *= z;
        result += coeffs[i] * power;
    }
    return result;
}

//------------------------------------------------------------------------------
// std::vector<std::complex<double>> durandKerner(const std::vector<std::complex<double>>& coeffs, int maxIter = 1000)
//------------------------------------------------------------------------------
/**
 * Finds the roots of a polynomial using the Durand-Kerner method.
 *
 * @param coeffs  Coefficients of the polynomial in descending order of degree.
 * @param maxIter Maximum number of iterations (default: 1000).
 * @return        Vector containing the complex roots of the polynomial.
 */
//------------------------------------------------------------------------------
std::vector<std::complex<double>> durandKerner(const std::vector<std::complex<double>>& coeffs, int maxIter = 1000) {
    size_t degree = coeffs.size() - 1;
    std::vector<std::complex<double>> roots(degree);

    // Initial guess for roots (using a circle in the complex plane)
    std::complex<double> initialGuess = std::polar(1.0, 2.0 * M_PI / degree);
    for (size_t i = 0; i < degree; ++i) {
        roots[i] = std::pow(initialGuess, static_cast<double>(i));
    }

    // Iterate to refine the root estimates
    for (int iter = 0; iter < maxIter; ++iter) {
        std::vector<std::complex<double>> newRoots = roots;
        bool converged = true;
        for (size_t i = 0; i < degree; ++i) {
            std::complex<double> product = 1.0;
            for (size_t j = 0; j < degree; ++j) {
                if (i != j) {
                    product *= (roots[i] - roots[j]);
                }
            }
            newRoots[i] = roots[i] - evaluatePolynomial(coeffs, roots[i]) / product;
            if (std::abs(newRoots[i] - roots[i]) > TOL) {
                converged = false;
            }
        }
        roots = newRoots;
        if (converged) {
            break;
        }
    }
    return roots;
}

//------------------------------------------------------------------------------
// bool isReal(const std::complex<double>& num)
//------------------------------------------------------------------------------
/**
 * Checks if a complex number is effectively real within a defined tolerance.
 *
 * @param num The complex number to be checked.
 * @return True if the imaginary part of the complex number is within the tolerance, indicating it is effectively real; otherwise, false.
 */
//------------------------------------------------------------------------------
bool isReal(const std::complex<double>& num) {
    return std::abs(num.imag()) < TOL;
}

//------------------------------------------------------------------------------
// void largestRoot(double c1, double c2, double c3, double c4, double c5, double c6, double c7, double c8, double c9)
//------------------------------------------------------------------------------
/**
 * Finds the largest real root of a polynomial of degree 8 using the Durand-Kerner method.
 *
 * @param c1-c9 Coefficients of the polynomial in descending order of degree.
 * @return The largest real root of the polynomial.
 */
//------------------------------------------------------------------------------
double largestRoot(double c1, double c2, double c3, double c4, double c5, double c6, double c7, double c8, double c9) {
    std::vector<std::complex<double>> coeffs = {
            c1, c2, c3, c4, c5, c6, c7, c8, c9
    };

    // Find the roots of the polynomial
    std::vector<std::complex<double>> rootarr = durandKerner(coeffs);

    // Initialize bigr2 to a very small number
    double bigr2 = -99999990.0;

    // Find the largest real root
    for (const auto& root : rootarr) {
        if (root.real() > bigr2 && isReal(root)) {
            bigr2 = root.real();
        }
    }

    return bigr2;
}