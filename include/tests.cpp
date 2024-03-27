#include <iostream>
#include <cassert>
#include "Matrix.cpp"

void test_creacion_matrix() {
    Matrix m(2, 2);

    assert(m.numFilas() == 2);
    assert(m.numColumnas() == 2);

    m.set(0, 0, 1.0);
    m.set(0, 1, 2.0);
    m.set(1, 0, 3.0);
    m.set(1, 1, 4.0);

    assert(m.get(0, 0) == 1.0);
    assert(m.get(0, 1) == 2.0);
    assert(m.get(1, 0) == 3.0);
    assert(m.get(1, 1) == 4.0);
}

void test_operaciones_aritmeticas() {
    Matrix m1(2, 2);
    m1.set(0, 0, 1.0);
    m1.set(0, 1, 2.0);
    m1.set(1, 0, 3.0);
    m1.set(1, 1, 4.0);

    Matrix m2(2, 2);
    m2.set(0, 0, 5.0);
    m2.set(0, 1, 6.0);
    m2.set(1, 0, 7.0);
    m2.set(1, 1, 8.0);

    // Suma
    Matrix suma = m1 + m2;
    assert(suma.get(0, 0) == 6.0);
    assert(suma.get(0, 1) == 8.0);
    assert(suma.get(1, 0) == 10.0);
    assert(suma.get(1, 1) == 12.0);

    // Resta
    Matrix resta = m2 - m1;
    assert(resta.get(0, 0) == 4.0);
    assert(resta.get(0, 1) == 4.0);
    assert(resta.get(1, 0) == 4.0);
    assert(resta.get(1, 1) == 4.0);

    // Multiplicación por escalar
    Matrix multiplicacion = m1 * 2.0;
    assert(multiplicacion.get(0, 0) == 2.0);
    assert(multiplicacion.get(0, 1) == 4.0);
    assert(multiplicacion.get(1, 0) == 6.0);
    assert(multiplicacion.get(1, 1) == 8.0);
}

void test_transposicion() {
    Matrix m(2, 3);
    m.set(0, 0, 1.0);
    m.set(0, 1, 2.0);
    m.set(0, 2, 3.0);
    m.set(1, 0, 4.0);
    m.set(1, 1, 5.0);
    m.set(1, 2, 6.0);

    Matrix transpuesta = m.transpose();
    assert(transpuesta.numFilas() == 3);
    assert(transpuesta.numColumnas() == 2);
    assert(transpuesta.get(0, 0) == 1.0);
    assert(transpuesta.get(0, 1) == 4.0);
    assert(transpuesta.get(1, 0) == 2.0);
    assert(transpuesta.get(1, 1) == 5.0);
    assert(transpuesta.get(2, 0) == 3.0);
    assert(transpuesta.get(2, 1) == 6.0);
}

void test_multiplicacion_matrices() {
    Matrix m1(2, 3);
    m1.set(0, 0, 1.0);
    m1.set(0, 1, 2.0);
    m1.set(0, 2, 3.0);
    m1.set(1, 0, 4.0);
    m1.set(1, 1, 5.0);
    m1.set(1, 2, 6.0);

    Matrix m2(3, 2);
    m2.set(0, 0, 7.0);
    m2.set(0, 1, 8.0);
    m2.set(1, 0, 9.0);
    m2.set(1, 1, 10.0);
    m2.set(2, 0, 11.0);
    m2.set(2, 1, 12.0);

    Matrix resultado = m1 * m2;
    assert(resultado.numFilas() == 2);
    assert(resultado.numColumnas() == 2);
    assert(resultado.get(0, 0) == 58.0);
    assert(resultado.get(0, 1) == 64.0);
    assert(resultado.get(1, 0) == 139.0);
    assert(resultado.get(1,1) == 154.0);
}

int main() {
    // Ejecutar los tests
    test_creacion_matrix();
    std::cout << "Test de creación de matriz pasado." << std::endl;

    test_operaciones_aritmeticas();
    std::cout << "Test de operaciones aritméticas pasado." << std::endl;

    test_transposicion();
    std::cout << "Test de transposición pasado." << std::endl;

    test_multiplicacion_matrices();
    std::cout << "Test de multiplicación de matrices pasado." << std::endl;

    return 0;
}

