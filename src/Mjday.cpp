#include "../include/Mjday.h"

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
double Mjday(int year, int month, int day, int hour, int min, double sec){
    double dyear = static_cast<double>(year);
    double dmonth = static_cast<double>(month);
    double dday = static_cast<double>(day);
    double dhour = static_cast<double>(hour);
    double dmin = static_cast<double>(min);

    double jd = 367.0 * year - floor( (7.0 * (year + floor( (month + 9.0) / 12.0) ) ) * 0.25 ) + floor( 275.0 * month / 9.0 ) + day + 1721013.5 + ( (sec/60.0 + min ) / 60.0 + hour ) / 24.0;

    return jd-2400000.5;
}

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
double Mjday(int year, int month, int day){
    return Mjday(year, month, day, 0, 0, 0);
}