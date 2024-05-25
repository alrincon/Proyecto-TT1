#include "../include/Legendre.h"
#include <stdio.h>
#include <vector>
#include <iostream>

//------------------------------------------------------------------------------
// void Legendre(int n, int m, double fi, Matrix& pnm, Matrix& dpnm)
//------------------------------------------------------------------------------
/**
 * Computes the associated Legendre functions and their derivatives up to degree n and order m at a given latitude.
 *
 * @param n     Maximum degree of the Legendre functions.
 * @param m     Maximum order of the Legendre functions.
 * @param fi    Latitude in radians at which the Legendre functions are evaluated.
 * @param pnm   Matrix storing the Legendre functions up to degree n and order m.
 * @param dpnm  Matrix storing the derivatives of the Legendre functions up to degree n and order m.
 */
//------------------------------------------------------------------------------
void Legendre(int n, int m, double fi, Matrix& pnm, Matrix& dpnm) {
    pnm(1, 1) = 1;
    dpnm(1, 1) = 0;
    pnm(2, 2) = sqrt(3.0) * cos(fi);
    dpnm(2, 2) = -sqrt(3.0) * sin(fi);

    //diagonal coefficients
    for (int i = 2; i <= n; i++) {
        pnm(i + 1, i + 1) = sqrt((2.0 * i + 1.0) / (2.0 * i)) * cos(fi) * pnm(i, i);
    }

    for (int i = 2; i <= n; i++) {
        //dpnm(i + 1, i + 1) = sqrt((2.0 * i + 1.0) / (2.0 * i)) * ((cos(fi) * dpnm(i, i)) - (sin(fi) * pnm(i, i)));
        dpnm(i+1,i+1) = sqrt((2.0*i+1.0)/(2.0*i))*((cos(fi)*dpnm(i,i)) - (sin(fi)*pnm(i,i)));
    }

    // horizontal first step coefficients
    for (int i = 1; i <= n; i++) {
        pnm(i + 1, i) = sqrt(2.0 * i + 1.0) * sin(fi) * pnm(i, i);
    }

    for (int i = 1; i <= n; i++) {
        //dpnm(i + 1, i) = sqrt(2.0 * i + 1.0) * ((cos(fi) * pnm(i, i)) + (sin(fi) * dpnm(i, i)));
        dpnm(i+1,i)= sqrt(2.0*i+1.0)*((cos(fi)*pnm(i,i))+(sin(fi)*dpnm(i,i)));

    }

    // horizontal second step coefficients
    int j = 0;
    int k = 2;

    while (1) {
        for(int i = k; i <= n; i++) {
            pnm(i + 1, j + 1) = sqrt((2.0 * i + 1.0) / (1.0*(i - j) * (i + j))) *((sqrt(2.0 * i - 1.0) * sin(fi) * pnm(i, j + 1)) - (sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0)) * pnm(i - 1, j + 1)));
        }

        j = j + 1;
        k = k + 1;

        if (j > m){
            break;
        }
    }

    j = 0;
    k = 2;

    while (1) {
        for(int i = k; i <= n; i++) {
            //dpnm(i+1,j+1)=sqrt((2.0*i+1.9)/(1.0*(i-j)*(i+j)))*((sqrt(2.0*i-1.0)*sin(fi)*dpnm(i,j+1))+(sqrt(2.0*i-1.0)*cos(fi)*pnm(i,j+1))-(sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0))*dpnm(i-1,j+1)));
            dpnm(i+1,j+1)=sqrt((2.0*i+1.0)/((i-j)*(i+j)))*((sqrt(2.0*i-1.0)*sin(fi)*dpnm(i,j+1))+(sqrt(2.0*i-1.0)*cos(fi)*pnm(i,j+1))-(sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0))*dpnm(i-1,j+1)));

        }

        j = j + 1;
        k = k + 1;

        if (j > m){
            break;
        }
    }
}