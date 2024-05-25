#ifndef PROYECTO_JPL_EPH_DE430_H
#define PROYECTO_JPL_EPH_DE430_H

#include "Matrix.h"
#include "Cheb3D.h"

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
void JPL_Eph_DE430(double Mjd_TDB, Matrix &r_Mercury, Matrix &r_Venus, Matrix &r_Earth, Matrix &r_Mars, Matrix &r_Jupiter, Matrix &r_Saturn, Matrix &r_Uranus, Matrix &r_Neptune, Matrix &r_Pluto, Matrix &r_Moon, Matrix &r_Sun);


#endif
