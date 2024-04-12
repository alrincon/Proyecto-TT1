#include "../include/Matrix.h"
#include <assert.h>
#include <iostream>
using namespace std;

void test_constructor_and_destructor() {
    // Test constructor básico
    Matrix m1(2, 3);
    assert(m1(0, 0) == 0 && "Constructor básico no inicializa a cero.");

    // Test constructor con valores
    double vals[] = {1, 2, 3, 4, 5, 6};
    Matrix m2(2, 3, vals, 6);
    assert(m2(0, 0) == 1 && m2(0, 1) == 2 && m2(1, 2) == 6 && "Constructor con valores no funciona correctamente.");

    // Test de destructor implícito por la falta de memory leaks o errores en tiempo de ejecución
    {
        Matrix m3(2, 2);
    }
}

void test_copy_constructor_and_assignment() {
    double vals[] = {1, 2, 3, 4};
    Matrix m1(2, 2, vals, 4);
    Matrix m2 = m1; // Test copy constructor
    assert(m2(0, 0) == 1 && m2(1, 1) == 4 && "Copy constructor no funciona correctamente.");

    Matrix m3(2, 2);
    m3 = m2; // Test operator=
    assert(m3(0, 0) == 1 && m3(1, 1) == 4 && "Assignment operator no funciona correctamente.");
}

void test_arithmetic_operations() {
    double vals1[] = {1, 2, 3, 4};
    double vals2[] = {4, 3, 2, 1};
    Matrix m1(2, 2, vals1, 4);
    Matrix m2(2, 2, vals2, 4);

    Matrix sum = m1 + m2;
    assert(sum(0, 0) == 5 && sum(1, 1) == 5 && "Addition operator no funciona correctamente.");

    Matrix diff = m1 - m2;
    assert(diff(0, 0) == -3 && diff(1, 1) == 3 && "Subtraction operator no funciona correctamente.");

    Matrix prod = m1 * m2;
    // Aquí suponemos un comportamiento específico del operador *, que debe ser validado según la implementación
    // por ejemplo, si m1 * m2 en (0,0) es 1*4 + 2*2 = 8
    assert(prod(0, 0) == 8 && "Multiplication operator no funciona correctamente.");
}

void test_bounds_checking() {
    Matrix m(2, 2);
    bool caught = false;
    try {
        m(2, 2);
    } catch (const std::out_of_range& e) {
        caught = true;
    }
    assert(caught && "Out-of-bounds access no atrapa excepciones como se esperaba.");
}

void test_all() {
    test_constructor_and_destructor();
    test_copy_constructor_and_assignment();
    test_arithmetic_operations();
    test_bounds_checking();
}
