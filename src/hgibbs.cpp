#include "../include/hgibbs.h"
#include "../include/angl.h"

//------------------------------------------------------------------------------
// void hgibbs (Matrix* r1, Matrix* r2, Matrix* r3,
//              double Mjd1, double Mjd2, double Mjd3,
//              Matrix& v2, double &theta, double &theta1,
//              double &copa, char* &error)
//------------------------------------------------------------------------------
/**
 * Implements the Herrick-Gibbs method for solving the three-body problem,
 * providing the velocity vector at the second position based on three
 * position vectors and their respective times.
 *
 * @param r1    Pointer to the position vector at the first time.
 * @param r2    Pointer to the position vector at the second time.
 * @param r3    Pointer to the position vector at the third time.
 * @param Mjd1  Modified Julian Date (MJD) of the first time.
 * @param Mjd2  Modified Julian Date (MJD) of the second time.
 * @param Mjd3  Modified Julian Date (MJD) of the third time.
 * @param v2    Output velocity vector at the second position.
 * @param theta Angle between r1 and r2.
 * @param theta1 Angle between r2 and r3.
 * @param copa  Angle between the normal to the plane of r2 and r3 and r1.
 * @param error Error message if any.
 */
//------------------------------------------------------------------------------
void hgibbs (Matrix* r1, Matrix* r2, Matrix* r3, double Mjd1, double Mjd2, double Mjd3, Matrix& v2, double &theta, double &theta1,double &copa, char* &error){
    error = static_cast<char *>(malloc(100));

    strcpy(error,"ok");
    theta = 0.0;
    theta1= 0.0;

    double magr1 = (*r1).norm();
    double magr2 = (*r2).norm();
    double magr3 = (*r3).norm();

    v2(1,1)= 0.0;
    v2(2,1)= 0.0;
    v2(3,1)= 0.0;

    double tolangle = 0.01745329251994;
    double dt21 = (Mjd2-Mjd1)*86400.0;
    double dt31 = (Mjd3-Mjd1)*86400.0;
    double dt32 = (Mjd3-Mjd2)*86400.0;

    Matrix p = crossProduct( r2,r3 );
    Matrix pn = unit(&p);
    Matrix r1n = unit( r1 );
    copa =  asin( dotProduct( &pn,&r1n ) );

    if (abs( dotProduct(&r1n,&pn) ) > 0.017452406 ) {
        strcpy(error, "not coplanar");
    }

    theta  = angl( r1,r2 );
    theta1 = angl( r2,r3 );

    if ( (theta > tolangle) | (theta1 > tolangle) ) {
        strcpy(error, "angl > 1ø");
    }

    double term1 = -dt32*( 1.0/(dt21*dt31) + GM_Earth/(12.0*magr1*magr1*magr1) );
    double term2= (dt32-dt21)*( 1.0/(dt21*dt32) + GM_Earth/(12.0*magr2*magr2*magr2) );
    double term3=  dt21*( 1.0/(dt32*dt31) + GM_Earth/(12.0*magr3*magr3*magr3) );

    v2 = (*r1)*term1 + (*r2)*term2 + (*r3)*term3;
}