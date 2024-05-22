#include "../include/Mjday.h"

double Mjday(int year, int month, int day, int hour, int min, double sec){
    double dyear = static_cast<double>(year);
    double dmonth = static_cast<double>(month);
    double dday = static_cast<double>(day);
    double dhour = static_cast<double>(hour);
    double dmin = static_cast<double>(min);

    double jd = 367.0 * dyear - floor( (7.0 * (dyear + floor( (dmonth + 9.0) / 12.0) ) ) * 0.25 ) + floor( 275.0 * dmonth / 9.0 ) + dday + 1721013.5 + ( (sec/60.0 + dmin ) / 60.0 + dhour ) / 24.0;

    return jd-2400000.5;
}

double Mjday(int year, int month, int day){
    return Mjday(year, month, day, 0, 0, 0);
}