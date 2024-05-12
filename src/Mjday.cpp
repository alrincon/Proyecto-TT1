#include "../include/Mjday.h"

double Mjday(int year, int month, int day, int hour, int min, double sec){
    double jd = 367.0 * year - floor( (7 * (year + floor( (month + 9) / 12.0) ) ) * 0.25 ) + floor( 275 * month / 9.0 ) + day + 1721013.5 + ( (sec/60.0 + min ) / 60.0 + hour ) / 24.0;

    return jd-2400000.5;
}

double Mjday(int year, int month, int day){
    return Mjday(year, month, day, 0, 0, 0);
}