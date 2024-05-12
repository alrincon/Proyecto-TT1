#include "../include/Position.h"

Matrix Position(double lon, double lat, double h){
    Matrix res(3,1);

    double e2 = f_Earth*(2-f_Earth);
    double CosLat = cos(lat);
    double SinLat = sin(lat);

    double N = R_Earth / sqrt(1-e2*SinLat*SinLat);

    res(1,1) = (N+h)*CosLat*cos(lon);
    res(2,1) = (N+h)*CosLat*sin(lon);
    res(3,1) = ((1.0-e2)*N+h)*SinLat;

    return res;
}
