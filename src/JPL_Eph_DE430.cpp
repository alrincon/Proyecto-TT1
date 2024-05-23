#include "../include/JPL_Eph_DE430.h"


extern Matrix PC;

int findPC(double JD){
    for(int i = 1; i <= PC.getFilas(); i++){
        if(PC(i,1) <= JD && JD <= PC(i,2)){
            return i;
        }
    }

    return 1;
}

Matrix extractRow(Matrix* input, int column){
    Matrix res(input->getFilas(), 1);
    for(int i = 1; i <= input->getFilas(); i++){
        res(i, 1) = (*input)(i, column);
    }

    return res;
}

Matrix extractColumn(Matrix* input, int row){
    Matrix res(1, input->getColumnas());
    for(int i = 1; i <= input->getColumnas(); i++){
        res( 1,i) = (*input)(row, i);
    }

    return res;
}

Matrix extractVector(Matrix* input, int start, int finish){
    int n = finish-start;
    Matrix res(1, n+1);

    for(int i = 0; i <= n; i++){
        res(1,i+1) = (*input)(1,start + i);
    }

    return res;
}

Matrix rowSequence(int start, int finish, int step){
    int n = (finish-start)/step;
    Matrix res(1, n+1);

    for(int i = 0; i <= n; i++){
        res(1,i+1) = start + i*step;
    }

    return res;
}

Matrix concatenateVector(Matrix *v1, Matrix *v2){
    Matrix res(1, v1->getColumnas() + v2->getColumnas());

    for(int i = 1; i <= v1->getColumnas(); i++){
        res(1,i) = (*v1)(1,i);
    }

    for(int i = 1; i <= v2->getColumnas(); i++){
        res(1,v1->getColumnas()+i) = (*v1)(1,i);
    }

    return res;
}



void JPL_Eph_DE430(double Mjd_TDB, Matrix &r_Mercury, Matrix &r_Venus, Matrix &r_Earth, Matrix &r_Mars, Matrix &r_Jupiter, Matrix &r_Saturn, Matrix &r_Uranus, Matrix &r_Neptune, Matrix &r_Pluto, Matrix &r_Moon, Matrix &r_Sun) {
    double JD = Mjd_TDB + 2400000.5;
    int i = findPC(JD);
    Matrix PCtemp = extractColumn(&PC, i);

    double t1 = PCtemp(1, 1) - 2400000.5; // MJD at start of interval
    double dt = Mjd_TDB - t1;

    Matrix temp = rowSequence(231, 270, 13);

    Matrix Cx_Earth = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
    Matrix Cy_Earth = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
    Matrix Cz_Earth = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

    for (int k = 1; k <= temp.getColumnas(); k++) {
        temp(1, k) = temp(1, k) + 39;
    }


    Matrix Cx = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
    Matrix Cy = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
    Matrix Cz = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

    Matrix aux1 = concatenateVector(&Cx_Earth, &Cx);
    Cx_Earth.redefine(&aux1);

    Matrix aux2 = concatenateVector(&Cy_Earth, &Cy);
    Cy_Earth.redefine(&aux2);

    Matrix aux3 = concatenateVector(&Cz_Earth, &Cz);
    Cz_Earth.redefine(&aux3);

    int j = 0;
    double Mjd0;

    if (0 <= dt && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    } else {
        if (16 < dt && dt <= 32) {
            j = 1;
            Mjd0 = t1 + 16 * j;
        }
    }

    Matrix Cxt = extractVector(&Cx_Earth, 13 * j + 1, 13 * j + 13);
    Matrix Cyt = extractVector(&Cy_Earth, 13 * j + 1, 13 * j + 13);
    Matrix Czt = extractVector(&Cz_Earth, 13 * j + 1, 13 * j + 13);

    r_Earth = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 16, &Cxt, &Cyt, &Czt).transpose() * 1e3;

    //Moon
    temp = rowSequence(441, 480, 13);
    Matrix Cx_Moon = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
    Matrix Cy_Moon = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
    Matrix Cz_Moon = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

    for (int r = 1; r <= 7; r++) {
        for (int k = 1; k <= temp.getColumnas(); k++) {
            temp(1, k) = temp(1, k) + 39;
        }

        Cx = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
        Cy = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
        Cz = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

        Matrix aux1 = concatenateVector(&Cx_Moon, &Cx);
        Cx_Moon.redefine(&aux1);

        Matrix aux2 = concatenateVector(&Cy_Moon, &Cy);
        Cy_Moon.redefine(&aux2);

        Matrix aux3 = concatenateVector(&Cz_Moon, &Cz);
        Cz_Moon.redefine(&aux3);
    }

    if (0 <= dt && dt <= 4) {
        j = 0;
        Mjd0 = t1;
    } else if (4 < dt && dt <= 8) {
        j = 1;
        Mjd0 = t1 + 4 * j;
    } else if (8 < dt && dt <= 12) {
        j = 2;
        Mjd0 = t1 + 4 * j;
    } else if (12 < dt && dt <= 16) {
        j = 3;
        Mjd0 = t1 + 4 * j;
    } else if (16 < dt && dt <= 20) {
        j = 4;
        Mjd0 = t1 + 4 * j;
    } else if (20 < dt && dt <= 24) {
        j = 5;
        Mjd0 = t1 + 4 * j;
    } else if (24 < dt && dt <= 28) {
        j = 6;
        Mjd0 = t1 + 4 * j;
    } else if (28 < dt && dt <= 32) {
        j = 7;
        Mjd0 = t1 + 4 * j;
    }

    Cxt = extractVector(&Cx_Moon, 13 * j + 1, 13 * j + 13);
    Cyt = extractVector(&Cx_Moon, 13 * j + 1, 13 * j + 13);
    Czt = extractVector(&Cx_Moon, 13 * j + 1, 13 * j + 13);

    r_Moon = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 4, &Cxt, &Cyt, &Czt).transpose() * 1e3;

    //Sun
    temp = rowSequence(753, 786, 11);
    Matrix Cx_Sun = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
    Matrix Cy_Sun = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
    Matrix Cz_Sun = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

    for (int k = 1; k <= temp.getColumnas(); k++) {
        temp(1, k) = temp(1, k) + 33;
    }

    Cx = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
    Cy = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
    Cz = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

    Matrix aux21 = concatenateVector(&Cx_Sun, &Cx);
    Cx_Sun.redefine(&aux21);

    Matrix aux22 = concatenateVector(&Cy_Sun, &Cy);
    Cy_Sun.redefine(&aux22);

    Matrix aux23 = concatenateVector(&Cz_Sun, &Cz);
    Cz_Sun.redefine(&aux23);

    if (0 <= dt && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    } else if (16 < dt && dt <= 32) {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }

    Cxt = extractVector(&Cx_Sun, 11 * j + 1, 11 * j + 11);
    Cyt = extractVector(&Cx_Sun, 11 * j + 1, 11 * j + 11);
    Czt = extractVector(&Cx_Sun, 11 * j + 1, 11 * j + 11);

    r_Sun = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 16, &Cxt, &Cyt, &Czt).transpose() * 1e3;

    //Mercury
    temp = rowSequence(3, 45, 14);

    Matrix Cx_Mercury = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
    Matrix Cy_Mercury = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
    Matrix Cz_Mercury = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

    for (int r = 1; r <= 3; r++) {
        for (int k = 1; k <= temp.getColumnas(); k++) {
            temp(1, k) = temp(1, k) + 42;
        }

        Cx = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
        Cy = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
        Cz = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

        Matrix aux21 = concatenateVector(&Cx_Mercury, &Cx);
        Cx_Mercury.redefine(&aux21);

        Matrix aux22 = concatenateVector(&Cy_Mercury, &Cy);
        Cy_Mercury.redefine(&aux22);

        Matrix aux23 = concatenateVector(&Cz_Mercury, &Cz);
        Cz_Mercury.redefine(&aux23);
    }

    if (0 <= dt && dt <= 8) {
        j = 0;
        Mjd0 = t1;
    } else if (8 < dt && dt <= 16) {
        j = 1;
        Mjd0 = t1 + 8 * j;
    } else if (16 < dt && dt <= 24) {
        j = 2;
        Mjd0 = t1 + 8 * j;
    } else if (24 < dt && dt <= 32) {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }

    Cxt = extractVector(&Cx_Mercury, 14 * j + 1, 14 * j + 14);
    Cyt = extractVector(&Cx_Mercury, 11 * j + 1, 11 * j + 11);
    Czt = extractVector(&Cx_Mercury, 11 * j + 1, 11 * j + 11);

    r_Mercury = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0 + 8, &Cxt, &Cyt, &Czt).transpose() * 1e3;

    //Venus
    temp = rowSequence(171, 201, 10);
    Matrix Cx_Venus = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
    Matrix Cy_Venus = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
    Matrix Cz_Venus = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

    for (int k = 1; k <= temp.getColumnas(); k++) {
        temp(1, k) = temp(1, k) + 30;
    }

    Cx = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
    Cy = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
    Cz = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

    Matrix aux31 = concatenateVector(&Cx_Venus, &Cx);
    Cx_Venus.redefine(&aux31);

    Matrix aux32 = concatenateVector(&Cy_Venus, &Cy);
    Cy_Venus.redefine(&aux32);

    Matrix aux33 = concatenateVector(&Cz_Venus, &Cz);
    Cz_Venus.redefine(&aux33);

    if (0 <= dt && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    } else if (16 < dt && dt <= 32) {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }

    Cxt = extractVector(&Cx_Venus, 10 * j + 1, 10 * j + 10);
    Cyt = extractVector(&Cy_Venus, 10 * j + 1, 10 * j + 10);
    Czt = extractVector(&Cz_Venus, 10 * j + 1, 10 * j + 10);

    r_Venus = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 16, &Cxt, &Cyt, &Czt).transpose() * 1e3;

    //Mars
    temp = rowSequence(309, 342, 1);
    Matrix Cx_Mars = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
    Matrix Cy_Mars = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
    Matrix Cz_Mars = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

    j = 0;
    Mjd0 = t1;

    Cxt = extractVector(&Cx_Mars, 11 * j + 1, 11 * j + 11);
    Cyt = extractVector(&Cy_Mars, 11 * j + 1, 11 * j + 11);
    Czt = extractVector(&Cz_Mars, 11 * j + 1, 11 * j + 11);

    r_Mars = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 32, &Cxt, &Cyt, &Czt).transpose() * 1e3;

    //Jupiter

    temp = rowSequence(342, 366, 8);
    Matrix Cx_Jupiter = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
    Matrix Cy_Jupiter = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
    Matrix Cz_Jupiter = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

    j = 0;
    Mjd0 = t1;

    Cxt = extractVector(&Cx_Jupiter, 8 * j + 1, 8 * j + 8);
    Cyt = extractVector(&Cy_Jupiter, 8 * j + 1, 8 * j + 8);
    Czt = extractVector(&Cz_Jupiter, 8 * j + 1, 8 * j + 8);

    r_Jupiter = Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0 + 32, &Cxt, &Cyt, &Czt).transpose() * 1e3;

    //Saturn

    temp = rowSequence(366, 387, 7);
    Matrix Cx_Saturn = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
    Matrix Cy_Saturn = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
    Matrix Cz_Saturn = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

    j = 0;
    Mjd0 = t1;

    Cxt = extractVector(&Cx_Saturn, 7 * j + 1, 7 * j + 7);
    Cyt = extractVector(&Cy_Saturn, 7 * j + 1, 7 * j + 7);
    Czt = extractVector(&Cz_Saturn, 7 * j + 1, 7 * j + 7);

    r_Saturn = Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0 + 32, &Cxt, &Cyt, &Czt).transpose() * 1e3;

    //Uranus

    temp = rowSequence(387, 405, 6);
    Matrix Cx_Uranus = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
    Matrix Cy_Uranus = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
    Matrix Cz_Uranus = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

    j = 0;
    Mjd0 = t1;

    Cxt = extractVector(&Cx_Uranus, 6 * j + 1, 6 * j + 6);
    Cyt = extractVector(&Cy_Uranus, 6 * j + 1, 6 * j + 6);
    Czt = extractVector(&Cz_Uranus, 6 * j + 1, 6 * j + 6);

    r_Uranus = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, &Cxt, &Cyt, &Czt).transpose() * 1e3;

    //Neptune

    temp = rowSequence(405, 423, 6);
    Matrix Cx_Neptune = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
    Matrix Cy_Neptune = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
    Matrix Cz_Neptune = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

    j = 0;
    Mjd0 = t1;

    Cxt = extractVector(&Cx_Neptune, 6 * j + 1, 6 * j + 6);
    Cyt = extractVector(&Cy_Neptune, 6 * j + 1, 6 * j + 6);
    Czt = extractVector(&Cz_Neptune, 6 * j + 1, 6 * j + 6);

    r_Neptune = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, &Cxt, &Cyt, &Czt).transpose() * 1e3;

    //Pluto

    temp = rowSequence(423, 441, 6);
    Matrix Cx_Pluto = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
    Matrix Cy_Pluto = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
    Matrix Cz_Pluto = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

    j = 0;
    Mjd0 = t1;

    Cxt = extractVector(&Cx_Pluto, 6 * j + 1, 6 * j + 6);
    Cyt = extractVector(&Cy_Pluto, 6 * j + 1, 6 * j + 6);
    Czt = extractVector(&Cz_Pluto, 6 * j + 1, 6 * j + 6);

    r_Pluto = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, &Cxt, &Cyt, &Czt).transpose() * 1e3;

    //Nutations
    temp = rowSequence(819, 839, 10);

    Matrix Cx_Nutations = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
    Matrix Cy_Nutations = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);

    for (int r = 1; r <= 3; r++) {
        for (int k = 1; k <= temp.getColumnas(); k++) {
            temp(1, k) = temp(1, k) + 20;
        }

        Matrix Cx = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
        Matrix Cy = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);

        Cx_Nutations = concatenateVector(&Cx_Nutations, &Cx);
        Cy_Nutations = concatenateVector(&Cy_Nutations, &Cy);
    }

    if (0 <= dt && dt <= 8) {
        j = 0;
        Mjd0 = t1;
    } else if (8 < dt && dt <= 16) {
        j = 1;
        Mjd0 = t1 + 8 * j;
    } else if (16 < dt && dt <= 24) {
        j = 2;
        Mjd0 = t1 + 8 * j;
    } else if (24 < dt && dt <= 32) {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }

    Cxt = extractVector(&Cx_Nutations, 10 * j + 1, 10 * j + 10);
    Cyt = extractVector(&Cy_Nutations, 10 * j + 1, 10 * j + 10);
    Czt = Matrix(10, 1);

    Matrix Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 8, &Cxt, &Cyt, &Czt).transpose();

    //Librations

    temp = rowSequence(899, 929, 10);
    Matrix Cx_Librations = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
    Matrix Cy_Librations = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
    Matrix Cz_Librations = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

    for (int r = 1; r <= 3; r++) {
        for (int k = 1; k <= temp.getColumnas(); k++) {
            temp(1, k) = temp(1, k) + 30;
        }

        Matrix Cx = extractVector(&PCtemp, temp(1, 1), temp(1, 2) - 1);
        Matrix Cy = extractVector(&PCtemp, temp(1, 2), temp(1, 3) - 1);
        Matrix Cz = extractVector(&PCtemp, temp(1, 3), temp(1, 4) - 1);

        Cx_Librations = concatenateVector(&Cx_Librations, &Cx);
        Cy_Librations = concatenateVector(&Cy_Librations, &Cy);
        Cz_Librations = concatenateVector(&Cz_Librations, &Cz);
    }

    if (0<=dt && dt<=8) {
        j = 0;
        Mjd0 = t1;
    }else if(8<dt && dt<=16) {
        j = 1;
        Mjd0 = t1 + 8 * j;
    }else if (16<dt && dt<=24) {
        j = 2;
        Mjd0 = t1 + 8 * j;
    }else if(24<dt && dt<=32) {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }

    Cxt = extractVector(&Cx_Librations, 10*j+1,10*j+10);
    Cyt = extractVector(&Cy_Librations, 10*j+1,10*j+10);
    Czt = extractVector(&Cz_Librations, 10*j+1,10*j+10);

    Matrix Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, &Cxt, &Cyt, &Czt).transpose();

    double EMRAT = 81.30056907419062; // DE430
    double EMRAT1 = 1/(1+EMRAT);
    r_Earth = r_Earth-r_Moon*EMRAT1;
    r_Mercury = r_Earth*(-1)+r_Mercury;
    r_Venus = r_Earth*(-1)+r_Venus;
    r_Mars = r_Earth*(-1)+r_Mars;
    r_Jupiter = r_Earth*(-1)+r_Jupiter;
    r_Saturn = r_Earth*(-1)+r_Saturn;
    r_Uranus = r_Earth*(-1)+r_Uranus;
    r_Neptune = r_Earth*(-1)+r_Neptune;
    r_Pluto = r_Earth*(-1)+r_Pluto;
    r_Sun = r_Earth*(-1)+r_Sun;
}