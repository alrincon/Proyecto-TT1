#include "../include/elements.h"

void elements (Matrix* y, double &p, double &a, double &e, double &i, double &Omega, double &omega, double &M){
    double pi2 = 2*M_PI;

    Matrix r(1,3);                                      // Position
    r(1,1) = (*y)(1,1);
    r(1,2) = (*y)(1,2);
    r(1,3) = (*y)(1,3);

    Matrix v(1,3);                                      // Velocity
    v(1,1) = (*y)(1,4);
    v(1,2) = (*y)(1,5);
    v(1,3) = (*y)(1,6);

    Matrix h = crossProduct(&r,&v);                     // Areal velocity
    double magh = h.norm();
    p = magh*magh/GM_Earth;
    double H = h.norm();

    Omega = atan2 ( h(1,1), -h(1,2) );                                                  // Long. ascend. node
    Omega = realmod(Omega,pi2);
    i = atan2 ( sqrt(h(1,1)*h(1,1)+h(1,2)*h(1,2)), h(1,3) );             // Inclination
    double u = atan2 ( r(1,3)*H, -r(1,1)*h(1,2)+r(1,2)*h(1,1) );         // Arg. of latitude

    double R  = r.norm();                                                   // Distance

    a = 1/(2/R-dotProduct(&v,&v)/GM_Earth);                          // Semi-major axis

    double eCosE = 1-R/a;                                                   // e*cos(E)
    double eSinE = dotProduct(&r,&v)/sqrt(GM_Earth*a);               // e*sin(E)

    double e2 = eCosE*eCosE +eSinE*eSinE;
    e  = sqrt(e2);                                                          // Eccentricity
    double E  = atan2(eSinE,eCosE);                                         // Eccentric anomaly

    M  = realmod(E-eSinE,pi2);                                                 // Mean anomaly

    double nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);                        // True anomaly

    omega = realmod(u-nu,pi2);// Arg. of perihelion
}