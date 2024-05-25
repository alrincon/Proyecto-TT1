#ifndef PROYECTO_MJDAY_H
#define PROYECTO_MJDAY_H

#include <cmath>

//------------------------------------------------------------------------------
// double Mjday(int year, int month, int day, int hour, int min, double sec)
//------------------------------------------------------------------------------
/**
 * Computes the Modified Julian Date (MJD) given the date and time components.
 *
 * @param year  Year.
 * @param month Month.
 * @param day   Day.
 * @param hour  Hour.
 * @param min   Minute.
 * @param sec   Second.
 * @return      Modified Julian Date (MJD).
 */
//------------------------------------------------------------------------------
double Mjday(int year, int month, int day, int hour, int min, double sec);

//------------------------------------------------------------------------------
// double Mjday(int year, int month, int day)
//------------------------------------------------------------------------------
/**
 * Computes the Modified Julian Date (MJD) given the date components (time components set to zero).
 *
 * @param year  Year.
 * @param month Month.
 * @param day   Day.
 * @return      Modified Julian Date (MJD).
 */
//------------------------------------------------------------------------------
double Mjday(int year, int month, int day);


#endif //PROYECTO_MJDAY_H
