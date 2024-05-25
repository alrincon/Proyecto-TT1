#include "../include/doubler.h"

//------------------------------------------------------------------------------
// doubler(double cc1, double cc2, double magrsite1, double magrsite2, double magr1in, double magr2in, Matrix *los1, Matrix *los2, Matrix *los3, Matrix *rsite1, Matrix *rsite2, Matrix *rsite3, double t1, double t3, char direct, Matrix &r2, Matrix &r3, double &f1, double &f2, double &q1, double &magr1, double &magr2, double &a, double &deltae32)
//------------------------------------------------------------------------------
/**
 * Computes the vectors r2 and r3, as well as other parameters, for a double r approach.
 *
 * This function computes the vectors r2 and r3 for a double r approach, as well as other
 * parameters such as f1, f2, q1, magr1, magr2, a, and deltae32. It accepts several input
 * parameters and modifies some of them by reference.
 *
 * @param cc1 Double r coefficient for the first site.
 * @param cc2 Double r coefficient for the second site.
 * @param magrsite1 Magnitude of the position vector of the first site.
 * @param magrsite2 Magnitude of the position vector of the second site.
 * @param magr1in Magnitude of the input position vector r1.
 * @param magr2in Magnitude of the input position vector r2.
 * @param los1 Pointer to the line of sight unit vector for the first site.
 * @param los2 Pointer to the line of sight unit vector for the second site.
 * @param los3 Pointer to the line of sight unit vector for the third site.
 * @param rsite1 Pointer to the position vector of the first site.
 * @param rsite2 Pointer to the position vector of the second site.
 * @param rsite3 Pointer to the position vector of the third site.
 * @param t1 Time at the first site.
 * @param t3 Time at the third site.
 * @param direct Direct or indirect double r approach ('y' or 'n').
 * @param r2 Computed position vector r2.
 * @param r3 Computed position vector r3.
 * @param f1 Resulting time difference f1.
 * @param f2 Resulting time difference f2.
 * @param q1 Resulting time difference q1.
 * @param magr1 Magnitude of the computed position vector r1.
 * @param magr2 Magnitude of the computed position vector r2.
 * @param a Semi-major axis of the computed orbit.
 * @param deltae32 Difference between eccentric anomaly of orbits 3 and 2.
 *
 * @note The function modifies the values of r2, r3, f1, f2, q1, magr1, magr2, a, and deltae32.
 */
//------------------------------------------------------------------------------
void doubler(double cc1, double cc2, double magrsite1, double magrsite2, double magr1in, double magr2in, Matrix *los1, Matrix *los2, Matrix *los3, Matrix *rsite1, Matrix *rsite2, Matrix *rsite3, double t1, double t3, char direct, Matrix &r2, Matrix &r3, double &f1, double &f2, double &q1, double &magr1, double &magr2, double &a, double &deltae32){
    double rho1 = (-cc1 + sqrt(pow(cc1,2)-4.0*(pow(magrsite1,2)-pow(magr1in,2)))) / 2.0;
    double rho2 = (-cc2 + sqrt(pow(cc2,2)-4.0*(pow(magrsite2,2)-pow(magr2in,2)))) / 2.0;

    Matrix r1 = (*los1)*rho1 + *rsite1;
    r2 = (*los2)*rho2 + *rsite2;

    magr1 = r1.norm();
    magr2 = r2.norm();

    Matrix w(r1.getFilas(), r1.getColumnas());

    if (direct == 'y') {
        w = crossProduct(&r1, &r2)*(1/(magr1 * magr2));
    }else {
        w = crossProduct(&r1, &r2)*(-1/(magr1 * magr2));
    }

    double rho3 =  -dotProduct(rsite3,&w)/dotProduct(los3,&w);
    r3 = (*los3)*rho3 + *rsite3;
    double magr3 = r3.norm();

    double cosdv21 = dotProduct(&r2,&r1)/(magr2*magr1);
    double sindv21 = (crossProduct(&r2,&r1)).norm()/(magr2*magr1);
    double dv21 = atan2(sindv21,cosdv21);

    double cosdv31 = dotProduct(&r3,&r1)/(magr3*magr1);
    double sindv31 = sqrt(1.0 - pow(cosdv31,2));
    double dv31 = atan2(sindv31,cosdv31);

    double cosdv32 = dotProduct(&r3,&r2)/(magr3*magr2);
    double sindv32 = (crossProduct(&r3,&r2)).norm()/(magr3*magr2);
    double dv32 = atan2(sindv32,cosdv32);

    double c1;
    double c3;
    double p;

    if (dv31 > M_PI) {
        c1 = (magr2 * sindv32) / (magr1 * sindv31);
        c3 = (magr2 * sindv21) / (magr3 * sindv31);
        p = (c1 * magr1 + c3 * magr3 - magr2) / (c1 + c3 - 1);
    }else {
        c1 = (magr1 * sindv31) / (magr2 * sindv32);
        c3 = (magr1 * sindv21) / (magr3 * sindv32);
        p = (c3 * magr3 - c1 * magr2 + magr1) / (-c1 + c3 + 1);
    }

    double ecosv1 = p/magr1-1.0;
    double ecosv2 = p/magr2-1.0;
    double ecosv3 = p/magr3-1.0;

    double esinv2;

    if (abs(dv21 - M_PI) > 10e-6) {
        esinv2 = (-cosdv21 * ecosv2 + ecosv1) / sindv21;
    }else {
        esinv2 = (cosdv32 * ecosv2 - ecosv3) / sindv31;
    }

    double e = sqrt(pow(ecosv2,2)+pow(esinv2,2));
    a = p/(1-pow(e,2));

    double n;

    double s;
    double c;

    double sinde32;
    double cosde32;

    double sinde21;
    double cosde21;
    double deltae21;

    double deltam32;
    double deltam12;

    double sindh32;
    double sindh21;
    double deltah32;
    double deltah21;

    if (e*e < 0.99) {
        n = sqrt(GM_Earth / pow(a,3));

        s = magr2 / p * sqrt(1.0 - pow(e,2)) * esinv2;
        c = magr2 / p * (pow(e,2) + ecosv2);

        sinde32 = magr3 / sqrt(a * p) * sindv32 - magr3 / p * (1.0 - cosdv32) * s;
        cosde32 = 1.0 - magr2 * magr3 / (a * p) * (1 - cosdv32);
        deltae32 = atan2(sinde32, cosde32);

        sinde21 = magr1 / sqrt(a * p) * sindv21 + magr1 / p * (1.0 - cosdv21) * s;
        cosde21 = 1.0 - magr2 * magr1 / (a * p) * (1.0 - cosdv21);
        deltae21 = atan2(sinde21, cosde21);

        deltam32 = deltae32 + 2.0 * s * pow((sin(deltae32 / 2.0)), 2) - c * sin(deltae32);
        deltam12 = -deltae21 + 2.0 * s * pow((sin(deltae21 / 2.0)), 2) + c * sin(deltae21);
    }else {
        n = sqrt(GM_Earth / -pow(a,3));

        s = magr2 / p * sqrt(pow(e,2)*1.0 - 1.0) * esinv2;
        c = magr2 / p * (pow(e,2)*1.0 + 1.0*ecosv2);

        sindh32 = magr3 / sqrt(-a * p) * sindv32 - magr3 / p * (1 - cosdv32) * s;
        sindh21 = magr1 / sqrt(-a * p) * sindv21 + magr1 / p * (1 - cosdv21) * s;

        deltah32 = log(sindh32 + sqrt(pow(sindh32,2) + 1.0));
        deltah21 = log(sindh21 + sqrt(pow(sindh21,2) + 1.0));

        deltam32 = -deltah32 + 2 * s * pow((sinh(deltah32 / 2.0)),2) + c * sinh(deltah32);
        deltam12 = deltah21 + 2 * s * pow((sinh(deltah21 / 2.0)),2) - c * sinh(deltah21);

        deltae32 = deltah32;
    }

    f1 = t1-deltam12/n;
    f2 = t3-deltam32/n;

    q1 = sqrt(pow(f1,2)+pow(f2,2));
}
