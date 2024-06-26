#include "../include/gibbs.h"

//------------------------------------------------------------------------------
// gibbs(Matrix *r1, Matrix *r2, Matrix *r3, Matrix &v2, double &theta,
//       double &theta1, double &copa, char* &error)
//------------------------------------------------------------------------------
/**
 * Computes the velocity vector at the second point using the Gibbs method
 * given the position vectors at three points in space.
 *
 * @param r1 Pointer to the position vector at the first point.
 * @param r2 Pointer to the position vector at the second point.
 * @param r3 Pointer to the position vector at the third point.
 * @param v2 Reference to the velocity vector at the second point (output).
 * @param theta Reference to the first angular change (output).
 * @param theta1 Reference to the second angular change (output).
 * @param copa Reference to the angle between the plane of the three points
 *             and the first position vector (output).
 * @param error Reference to the error message (output).
 */
//------------------------------------------------------------------------------
void gibbs(Matrix *r1, Matrix *r2, Matrix *r3, Matrix &v2, double &theta, double &theta1, double &copa, char* &error){
    double small= 0.00000001;
    theta= 0.0;

    error = static_cast<char *>(malloc(100));

    strcpy(error,"ok");
    theta1= 0.0;

    double magr1 = (*r1).norm();
    double magr2 = (*r2).norm();
    double magr3 = (*r3).norm();

    v2(1,1)= 0.0;
    v2(2,1)= 0.0;
    v2(3,1)= 0.0;

    Matrix p(3,1);
    p = crossProduct( r2,r3 );
    Matrix q(3,1);
    q = crossProduct( r3,r1 );
    Matrix w(3,1);
    w = crossProduct( r1,r2 );

    Matrix pn(3, 1);
    pn = unit( &p );
    Matrix r1n(3, 1);

    r1n = unit( r1 );

    copa = asin( dotProduct(&pn,&r1n));

    if (abs(dotProduct(&r1n,&pn) ) > 0.017452406 ) {
        strcpy(error, "not coplanar");
    }

    Matrix d(3,1);
    d = p + q + w;
    double magd = d.norm();
    Matrix n(3,1);
    n = p*magr1 + q*magr2 + w*magr3;
    double magn = n.norm();
    Matrix nn(3,1);
    nn = unit( &n );
    Matrix dn(3,1);
    dn = unit( &d );

    // -------------------------------------------------------------
    //determine if  the orbit is possible. both d and n must be in
    //the same direction, and non-zero.
    //-------------------------------------------------------------

    if ( ( abs(magd)<small ) || ( abs(magn)<small ) || ( dotProduct(&nn,&dn) < small ) ) {
        strcpy(error, "impossible");
    }else {
        theta = angl(r1, r2);
        theta1 = angl(r2, r3);

        //-- -- -- -- -- -perform gibbs method to find v2-- -- -- -- -- -
        double r1mr2 = magr1 - magr2;
        double r3mr1 = magr3 - magr1;
        double r2mr3 = magr2 - magr3;
        Matrix s(3,1);
        s =  (*r3)*r1mr2 +  (*r2)*r3mr1 +  (*r1)*r2mr3;
        Matrix b(3,1);
        b = crossProduct(&d, r2);
        double l = sqrt(GM_Earth / (magd * magn));
        double tover2 = l / magr2;
        v2 =  b*tover2 +  s*l;
    }
}