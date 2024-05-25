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



//------------------------------------------------------------------------------
// findPC(double JD)
//------------------------------------------------------------------------------
/**
 * Calcula las posiciones de varios cuerpos celestes utilizando el efemérides DE430 de la JPL.
 *
 * @param Mjd_TDB   Fecha juliana modificada en el tiempo dinámico bariocéntrico (TDB).
 * @param r_Mercury Matriz para almacenar la posición de Mercurio.
 * @param r_Venus   Matriz para almacenar la posición de Venus.
 * @param r_Earth   Matriz para almacenar la posición de la Tierra.
 * @param r_Mars    Matriz para almacenar la posición de Marte.
 * @param r_Jupiter Matriz para almacenar la posición de Júpiter.
 * @param r_Saturn  Matriz para almacenar la posición de Saturno.
 * @param r_Uranus  Matriz para almacenar la posición de Urano.
 * @param r_Neptune Matriz para almacenar la posición de Neptuno.
 * @param r_Pluto   Matriz para almacenar la posición de Plutón.
 * @param r_Moon    Matriz para almacenar la posición de la Luna.
 * @param r_Sun     Matriz para almacenar la posición del Sol.
 */
//------------------------------------------------------------------------------
void JPL_Eph_DE430(double Mjd_TDB, Matrix &r_Mercury, Matrix &r_Venus, Matrix &r_Earth, Matrix &r_Mars, Matrix &r_Jupiter, Matrix &r_Saturn, Matrix &r_Uranus, Matrix &r_Neptune, Matrix &r_Pluto, Matrix &r_Moon, Matrix &r_Sun) {
    //solucionar uso compartido de matrices ¿funciones?

    double JD = Mjd_TDB + 2400000.5;
    int i = findPC(JD);
    Matrix PCtemp = extractColumn(&PC, i);

    double t1 = PCtemp(1, 1) - 2400000.5; // MJD at start of interval
    double dt = Mjd_TDB - t1;

    Matrix tempEarth = rowSequence(231, 270, 13);

    Matrix Cx_Earth = extractVector(&PCtemp, tempEarth(1, 1), tempEarth(1, 2) - 1);
    Matrix Cy_Earth = extractVector(&PCtemp, tempEarth(1, 2), tempEarth(1, 3) - 1);
    Matrix Cz_Earth = extractVector(&PCtemp, tempEarth(1, 3), tempEarth(1, 4) - 1);

    for (int k = 1; k <= tempEarth.getColumnas(); k++) {
        tempEarth(1, k) = tempEarth(1, k) + 39;
    }


    Matrix Cx = extractVector(&PCtemp, tempEarth(1, 1), tempEarth(1, 2) - 1);
    Matrix Cy = extractVector(&PCtemp, tempEarth(1, 2), tempEarth(1, 3) - 1);
    Matrix Cz = extractVector(&PCtemp, tempEarth(1, 3), tempEarth(1, 4) - 1);

    Matrix aux1 = concatenateVector(&Cx_Earth, &Cx);
    Cx_Earth.redefine(&aux1);

    Matrix aux2 = concatenateVector(&Cy_Earth, &Cy);
    Cy_Earth.redefine(&aux2);

    Matrix aux3 = concatenateVector(&Cz_Earth, &Cz);
    Cz_Earth.redefine(&aux3);

    int j = 0;
    double Mjd0 = 0;

    if (0 <= dt && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    } else {
        if (16 < dt && dt <= 32) {
            j = 1;
            Mjd0 = t1 + 16.0 * j;
        }
    }

    Matrix Cxt = extractVector(&Cx_Earth, 13 * j + 1, 13 * j + 13);
    Matrix Cyt = extractVector(&Cy_Earth, 13 * j + 1, 13 * j + 13);
    Matrix Czt = extractVector(&Cz_Earth, 13 * j + 1, 13 * j + 13);

    r_Earth = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 16, &Cxt, &Cyt, &Czt) * 1e3;

    //Moon
    Matrix tempMoon = rowSequence(441, 480, 13);
    Matrix Cx_Moon = extractVector(&PCtemp, tempMoon(1, 1), tempMoon(1, 2) - 1);
    Matrix Cy_Moon = extractVector(&PCtemp, tempMoon(1, 2), tempMoon(1, 3) - 1);
    Matrix Cz_Moon = extractVector(&PCtemp, tempMoon(1, 3), tempMoon(1, 4) - 1);

    for (int r = 1; r <= 7; r++) {
        for (int k = 1; k <= tempMoon.getColumnas(); k++) {
            tempMoon(1, k) = tempMoon(1, k) + 39;
        }

        Cx = extractVector(&PCtemp, tempMoon(1, 1), tempMoon(1, 2) - 1);
        Cy = extractVector(&PCtemp, tempMoon(1, 2), tempMoon(1, 3) - 1);
        Cz = extractVector(&PCtemp, tempMoon(1, 3), tempMoon(1, 4) - 1);

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

    Cxt.setTam(1, 13);
    Cyt.setTam(1, 13);
    Czt.setTam(1, 13);

    Cxt = extractVector(&Cx_Moon, 13 * j + 1, 13 * j + 13);
    Cyt = extractVector(&Cy_Moon, 13 * j + 1, 13 * j + 13);
    Czt = extractVector(&Cz_Moon, 13 * j + 1, 13 * j + 13);

    r_Moon = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 4, &Cxt, &Cyt, &Czt) * 1e3;

    //Sun
    Matrix tempSun = rowSequence(753, 786, 11);
    Matrix Cx_Sun = extractVector(&PCtemp, tempSun(1, 1), tempSun(1, 2) - 1);
    Matrix Cy_Sun = extractVector(&PCtemp, tempSun(1, 2), tempSun(1, 3) - 1);
    Matrix Cz_Sun = extractVector(&PCtemp, tempSun(1, 3), tempSun(1, 4) - 1);

    for (int k = 1; k <= tempSun.getColumnas(); k++) {
        tempSun(1, k) = tempSun(1, k) + 33;
    }

    Cx.setTam(1, tempSun(1, 2) - tempSun(1, 1));
    Cy.setTam(1, tempSun(1, 3) - tempSun(1, 2));
    Cz.setTam(1, tempSun(1, 4) - tempSun(1, 3));

    Cx = extractVector(&PCtemp, tempSun(1, 1), tempSun(1, 2) - 1);
    Cy = extractVector(&PCtemp, tempSun(1, 2), tempSun(1, 3) - 1);
    Cz = extractVector(&PCtemp, tempSun(1, 3), tempSun(1, 4) - 1);

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
        Mjd0 = t1 + 16.0 * j;
    }

    Cxt.setTam(1, 11);
    Cyt.setTam(1, 11);
    Czt.setTam(1, 11);
    Cxt = extractVector(&Cx_Sun, 11 * j + 1, 11 * j + 11);
    Cyt = extractVector(&Cy_Sun, 11 * j + 1, 11 * j + 11);
    Czt = extractVector(&Cz_Sun, 11 * j + 1, 11 * j + 11);

    r_Sun = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 16, &Cxt, &Cyt, &Czt) * 1e3;

    Matrix tempMercury = rowSequence(3, 45, 14);

    Matrix Cx_Mercury = extractVector(&PCtemp, tempMercury(1, 1), tempMercury(1, 2) - 1);
    Matrix Cy_Mercury = extractVector(&PCtemp, tempMercury(1, 2), tempMercury(1, 3) - 1);
    Matrix Cz_Mercury = extractVector(&PCtemp, tempMercury(1, 3), tempMercury(1, 4) - 1);


    for (int r = 1; r <= 3; r++) {
        for (int k = 1; k <= tempMercury.getColumnas(); k++) {
            tempMercury(1, k) = tempMercury(1, k) + 42;
        }

        Cx.setTam(1, tempMercury(1, 2) - tempMercury(1, 1));
        Cy.setTam(1, tempMercury(1, 3) - tempMercury(1, 2));
        Cz.setTam(1, tempMercury(1, 4) - tempMercury(1, 3));

        Cx = extractVector(&PCtemp, tempMercury(1, 1), tempMercury(1, 2) - 1);
        Cy = extractVector(&PCtemp, tempMercury(1, 2), tempMercury(1, 3) - 1);
        Cz = extractVector(&PCtemp, tempMercury(1, 3), tempMercury(1, 4) - 1);

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

    Cxt.setTam(1, 14);
    Cyt.setTam(1, 14);
    Czt.setTam(1, 14);

    Cxt = extractVector(&Cx_Mercury, 14 * j + 1, 14 * j + 14);
    Cyt = extractVector(&Cy_Mercury, 14 * j + 1, 14 * j + 14);
    Czt = extractVector(&Cz_Mercury, 14 * j + 1, 14 * j + 14);

    r_Mercury = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0 + 8, &Cxt, &Cyt, &Czt) * 1e3;

    //Venus
    Matrix tempVenus = rowSequence(171, 201, 10);
    Matrix Cx_Venus = extractVector(&PCtemp, tempVenus(1, 1), tempVenus(1, 2) - 1);
    Matrix Cy_Venus = extractVector(&PCtemp, tempVenus(1, 2), tempVenus(1, 3) - 1);
    Matrix Cz_Venus = extractVector(&PCtemp, tempVenus(1, 3), tempVenus(1, 4) - 1);

    for (int k = 1; k <= tempVenus.getColumnas(); k++) {
        tempVenus(1, k) = tempVenus(1, k) + 30;
    }

    Cx.setTam(1, tempVenus(1, 2) - tempVenus(1, 1));
    Cy.setTam(1, tempVenus(1, 3) - tempVenus(1, 2));
    Cz.setTam(1, tempVenus(1, 4) - tempVenus(1, 3));

    Cx = extractVector(&PCtemp, tempVenus(1, 1), tempVenus(1, 2) - 1);
    Cy = extractVector(&PCtemp, tempVenus(1, 2), tempVenus(1, 3) - 1);
    Cz = extractVector(&PCtemp, tempVenus(1, 3), tempVenus(1, 4) - 1);

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

    Cxt.setTam(1, 10);
    Cyt.setTam(1, 10);
    Czt.setTam(1, 10);

    Cxt = extractVector(&Cx_Venus, 10 * j + 1, 10 * j + 10);
    Cyt = extractVector(&Cy_Venus, 10 * j + 1, 10 * j + 10);
    Czt = extractVector(&Cz_Venus, 10 * j + 1, 10 * j + 10);

    r_Venus = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 16, &Cxt, &Cyt, &Czt) * 1e3;

    //Mars
    Matrix tempMars = rowSequence(309, 342, 11);
    Matrix Cx_Mars = extractVector(&PCtemp, tempMars(1, 1), tempMars(1, 2) - 1);
    Matrix Cy_Mars = extractVector(&PCtemp, tempMars(1, 2), tempMars(1, 3) - 1);
    Matrix Cz_Mars = extractVector(&PCtemp, tempMars(1, 3), tempMars(1, 4) - 1);

    j = 0;
    Mjd0 = t1;


    Cxt.setTam(1, 11);
    Cyt.setTam(1, 11);
    Czt.setTam(1, 11);

    Cxt = extractVector(&Cx_Mars, 11 * j + 1, 11 * j + 11);
    Cyt = extractVector(&Cy_Mars, 11 * j + 1, 11 * j + 11);
    Czt = extractVector(&Cz_Mars, 11 * j + 1, 11 * j + 11);

    r_Mars = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 32, &Cxt, &Cyt, &Czt) * 1e3;

    //Jupiter

    Matrix tempJupiter = rowSequence(342, 366, 8);
    Matrix Cx_Jupiter = extractVector(&PCtemp, tempJupiter(1, 1), tempJupiter(1, 2) - 1);
    Matrix Cy_Jupiter = extractVector(&PCtemp, tempJupiter(1, 2), tempJupiter(1, 3) - 1);
    Matrix Cz_Jupiter = extractVector(&PCtemp, tempJupiter(1, 3), tempJupiter(1, 4) - 1);

    j = 0;
    Mjd0 = t1;

    Cxt.setTam(1, 8);
    Cyt.setTam(1, 8);
    Czt.setTam(1, 8);

    Cxt = extractVector(&Cx_Jupiter, 8 * j + 1, 8 * j + 8);
    Cyt = extractVector(&Cy_Jupiter, 8 * j + 1, 8 * j + 8);
    Czt = extractVector(&Cz_Jupiter, 8 * j + 1, 8 * j + 8);

    r_Jupiter = Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0 + 32, &Cxt, &Cyt, &Czt) * 1e3;

    //Saturn

    Matrix tempSaturn = rowSequence(366, 387, 7);
    Matrix Cx_Saturn = extractVector(&PCtemp, tempSaturn(1, 1), tempSaturn(1, 2) - 1);
    Matrix Cy_Saturn = extractVector(&PCtemp, tempSaturn(1, 2), tempSaturn(1, 3) - 1);
    Matrix Cz_Saturn = extractVector(&PCtemp, tempSaturn(1, 3), tempSaturn(1, 4) - 1);

    j = 0;
    Mjd0 = t1;

    Cxt.setTam(1, 7);
    Cyt.setTam(1, 7);
    Czt.setTam(1, 7);

    Cxt = extractVector(&Cx_Saturn, 7 * j + 1, 7 * j + 7);
    Cyt = extractVector(&Cy_Saturn, 7 * j + 1, 7 * j + 7);
    Czt = extractVector(&Cz_Saturn, 7 * j + 1, 7 * j + 7);

    r_Saturn = Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0 + 32, &Cxt, &Cyt, &Czt) * 1e3;

    //Uranus

    Matrix tempUranus = rowSequence(387, 405, 6);
    Matrix Cx_Uranus = extractVector(&PCtemp, tempUranus(1, 1), tempUranus(1, 2) - 1);
    Matrix Cy_Uranus = extractVector(&PCtemp, tempUranus(1, 2), tempUranus(1, 3) - 1);
    Matrix Cz_Uranus = extractVector(&PCtemp, tempUranus(1, 3), tempUranus(1, 4) - 1);

    j = 0;
    Mjd0 = t1;

    Cxt.setTam(1, 6);
    Cyt.setTam(1, 6);
    Czt.setTam(1, 6);

    Cxt = extractVector(&Cx_Uranus, 6 * j + 1, 6 * j + 6);
    Cyt = extractVector(&Cy_Uranus, 6 * j + 1, 6 * j + 6);
    Czt = extractVector(&Cz_Uranus, 6 * j + 1, 6 * j + 6);

    r_Uranus = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, &Cxt, &Cyt, &Czt) * 1e3;

    //Neptune

    Matrix tempNeptune = rowSequence(405, 423, 6);
    Matrix Cx_Neptune = extractVector(&PCtemp, tempNeptune(1, 1), tempNeptune(1, 2) - 1);
    Matrix Cy_Neptune = extractVector(&PCtemp, tempNeptune(1, 2), tempNeptune(1, 3) - 1);
    Matrix Cz_Neptune = extractVector(&PCtemp, tempNeptune(1, 3), tempNeptune(1, 4) - 1);

    j = 0;
    Mjd0 = t1;

    Cxt = extractVector(&Cx_Neptune, 6 * j + 1, 6 * j + 6);
    Cyt = extractVector(&Cy_Neptune, 6 * j + 1, 6 * j + 6);
    Czt = extractVector(&Cz_Neptune, 6 * j + 1, 6 * j + 6);

    r_Neptune = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, &Cxt, &Cyt, &Czt) * 1e3;

    //Pluto

    Matrix tempPluto = rowSequence(423, 441, 6);
    Matrix Cx_Pluto = extractVector(&PCtemp, tempPluto(1, 1), tempPluto(1, 2) - 1);
    Matrix Cy_Pluto = extractVector(&PCtemp, tempPluto(1, 2), tempPluto(1, 3) - 1);
    Matrix Cz_Pluto = extractVector(&PCtemp, tempPluto(1, 3), tempPluto(1, 4) - 1);

    j = 0;
    Mjd0 = t1;

    Cxt = extractVector(&Cx_Pluto, 6 * j + 1, 6 * j + 6);
    Cyt = extractVector(&Cy_Pluto, 6 * j + 1, 6 * j + 6);
    Czt = extractVector(&Cz_Pluto, 6 * j + 1, 6 * j + 6);

    r_Pluto = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, &Cxt, &Cyt, &Czt) * 1e3;

    //Nutations
    Matrix tempNutations = rowSequence(819, 839, 10);

    Matrix Cx_Nutations = extractVector(&PCtemp, tempNutations(1, 1), tempNutations(1, 2) - 1);
    Matrix Cy_Nutations = extractVector(&PCtemp, tempNutations(1, 2), tempNutations(1, 3) - 1);

    for (int r = 1; r <= 3; r++) {
        for (int k = 1; k <= tempNutations.getColumnas(); k++) {
            tempNutations(1, k) = tempNutations(1, k) + 20;
        }

        Cx.setTam(1, tempNutations(1, 2) - tempNutations(1, 1));
        Cy.setTam(1, tempNutations(1, 3) - tempNutations(1, 2));

        Cx = extractVector(&PCtemp, tempNutations(1, 1), tempNutations(1, 2) - 1);
        Cy = extractVector(&PCtemp, tempNutations(1, 2), tempNutations(1, 3) - 1);

        Matrix aux1 = concatenateVector(&Cx_Nutations, &Cx);
        Cx_Nutations.redefine(&aux1);

        Matrix aux2 = concatenateVector(&Cy_Nutations, &Cy);
        Cy_Nutations.redefine(&aux2);

        /*Cx_Nutations = concatenateVector(&Cx_Nutations, &Cx);
        Cy_Nutations = concatenateVector(&Cy_Nutations, &Cy);*/
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

    Cxt.setTam(1, 10);
    Cyt.setTam(1, 10);
    Czt.setTam(1, 10);

    Cxt = extractVector(&Cx_Nutations, 10 * j + 1, 10 * j + 10);
    Cyt = extractVector(&Cy_Nutations, 10 * j + 1, 10 * j + 10);
    Czt = Matrix(1, 10);

    Matrix Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 8, &Cxt, &Cyt, &Czt);

    //Librations

    Matrix tempLibrations = rowSequence(899, 929, 10);
    Matrix Cx_Librations = extractVector(&PCtemp, tempLibrations(1, 1), tempLibrations(1, 2) - 1);
    Matrix Cy_Librations = extractVector(&PCtemp, tempLibrations(1, 2), tempLibrations(1, 3) - 1);
    Matrix Cz_Librations = extractVector(&PCtemp, tempLibrations(1, 3), tempLibrations(1, 4) - 1);

    for (int r = 1; r <= 3; r++) {
        for (int k = 1; k <= tempLibrations.getColumnas(); k++) {
            tempLibrations(1, k) = tempLibrations(1, k) + 30;
        }

        Matrix Cx = extractVector(&PCtemp, tempLibrations(1, 1), tempLibrations(1, 2) - 1);
        Matrix Cy = extractVector(&PCtemp, tempLibrations(1, 2), tempLibrations(1, 3) - 1);
        Matrix Cz = extractVector(&PCtemp, tempLibrations(1, 3), tempLibrations(1, 4) - 1);

        /*Cx_Librations = concatenateVector(&Cx_Librations, &Cx);
        Cy_Librations = concatenateVector(&Cy_Librations, &Cy);
        Cz_Librations = concatenateVector(&Cz_Librations, &Cz);*/

        Matrix aux1 = concatenateVector(&Cx_Librations, &Cx);
        Cx_Librations.redefine(&aux1);

        Matrix aux2 = concatenateVector(&Cy_Librations, &Cy);
        Cy_Librations.redefine(&aux2);

        Matrix aux3 = concatenateVector(&Cz_Librations, &Cz);
        Cz_Librations.redefine(&aux2);
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

    Matrix Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, &Cxt, &Cyt, &Czt);

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