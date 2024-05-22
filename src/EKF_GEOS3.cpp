#include "../include/Global.h"
#include "../include/SAT_const.h"
#include <stdio.h>
#include <string.h>
#include "../include/Mjday.h"
#include "../include/Position.h"
#include "../include/anglesg.h"
#include "../include/doubler.h"
#include <iostream>
#include <fstream>
#include "../include/DEInteg.h"
#include "../include/Accel.h"
#include "../include/VarEqn.h"
#include "../include/TimeUpdate.h"
#include "../include/AzElPa.h"
#include "../include/MeasUpdate.h"


#include "../include/EFK_Tests.h"


extern Matrix eopdata;
extern Matrix obs;
extern Matrix Cnm;
extern Matrix Snm;
extern aux AuxParam;


int main(){
    int nobs = 46;
    int infFile = 100;

    //fichero DE430Coeff.txt
    FILE *f1 = fopen("../data/DE430Coeff.txt", "r");
    if (!f1) {
        std::cerr << "No se pudo abrir el archivo." << std::endl;
        return 1;
    }

    for (int i = 1; i <= 2285; ++i) {
        for (int j = 1; j <= 1020; ++j) {
            if (fscanf(f1, "%lf", &PC(i, j)) != 1) {
                std::cerr << "Error al leer el archivo." << std::endl;
                fclose(f1);
                return 1;
            }
        }
    }

    fclose(f1);

    //fichero GGM03S.txt
    FILE *fp2 = fopen("../data/GGM03S.txt", "r");
    if (fp2 == nullptr) {
        printf("error");
        exit(1);
    }

    for (int i = 0; i < 181; i++) {
        int valor1;
        int valor2;
        double valor5;
        double valor6;

        for(int j = 0; j <= i; j++){
            fscanf(fp2, "%d %d %lf %lf %lf %lf", &valor1, &valor2, &Cnm(i+1, j+1), &Snm(i+1, j+1), &valor5, &valor6);
        }
    }

    if (fclose(fp2))   // Close the stream.
        perror("fclose error");
    else printf("File mylib/myfile closed successfully.\n");

    //parametros de orientación de la tierra
    FILE *fp3 = fopen("../data/eop19620101.txt", "r");
    if (fp3 == nullptr) {
        printf("error");
        exit(1);
    }

    for (int i = 1; i <= infFile; i++) {
        fscanf(fp3, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &eopdata(i, 1), &eopdata(i, 2), &eopdata(i, 3),
               &eopdata(i, 4), &eopdata(i, 5), &eopdata(i, 6),
               &eopdata(i, 7), &eopdata(i, 8), &eopdata(i, 9),
               &eopdata(i, 10), &eopdata(i, 11), &eopdata(i, 12),
               &eopdata(i, 13));
    }
    if (fclose(fp3))   // Close the stream.
        perror("fclose error");
    else printf("File mylib/myfile closed successfully.\n");


    //observaciones
    FILE *fp4 = fopen("../data/GEOS3.txt", "r");
    if (fp4 == nullptr) {
        printf("error");
        exit(1);
    }
    int year;
    int month;
    int day;

    int hour;
    int min;
    double sec;

    double az;
    double el;
    double Dist;

    int scancount = 1;
    while (fscanf(fp4, "%d/%d/%d %d:%d:%lf %lf %lf %lf",&year, &month, &day, &hour, &min, &sec,&az, &el, &Dist) == 9 && scancount <= nobs) {
        obs(scancount,1) = Mjday(year, month, day, hour,min, sec);
        obs(scancount,2) = Rad*az;
        obs(scancount,3) = Rad*el;
        obs(scancount,4) = 1e3*Dist;
        scancount++;
    }

    if (fclose(fp4))   // Close the stream.
        perror("fclose error");
    else printf("File mylib/myfile closed successfully.\n");

    double sigma_range = 92.5;          // [m]
    double sigma_az = 0.0224*Rad; // [rad]
    double sigma_el = 0.0139*Rad; // [rad]

    // Kaena Point station
    double lat = Rad*21.5748;     // [rad]
    double lon = Rad*(-158.2706); // [rad]
    double alt = 300.20;                // [m]

    Matrix Rs(3,1);
    Rs = Position(lon, lat, alt);


    double Mjd1 = obs(1,1);
    double Mjd2 = obs(9,1);
    double Mjd3 = obs(18,1);

    Matrix r2(1,3);
    Matrix v2(1,3);

    anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),Mjd1,Mjd2,Mjd3,&Rs,&Rs,&Rs, r2, v2);

    Matrix Y0_apr(1,6); //OJO

    Y0_apr(1,1) = r2(1,1);
    Y0_apr(1,2) = r2(1,2);
    Y0_apr(1,3) = r2(1,3);
    Y0_apr(1,4) = v2(1,1);
    Y0_apr(1,5) = v2(1,2);
    Y0_apr(1,6) = v2(1,3);


    double Mjd0 = Mjday(1995,1,29,02,38,0);

    double Mjd_UTC = obs(9,1);

    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;
    //parametros de modelo

    int n_eqn  = 6;

    Matrix Y(1,6);
    Matrix Y_old(1,6);

    cout << "llego" << endl;
    cout << Mjd0 << endl;
    cout << -(obs(9,1)-Mjd0)*86400.0 << endl;

    Y = DEInteg(&Accel,0.0,-(obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,&Y0_apr);

    cout << "llego" << endl;
    Matrix P(6,6);

    for (int i = 1; i <= 3; i++) {
        P(i, i) = 1e8;
    }

    for (int i = 4; i <= 6; i++) {
        P(i, i) = 1e3;
    }

    Matrix LT(3,3);
    LT = LTC(lon,lat);

    Matrix yPhi(42,1);
    Matrix Phi(6,6);

    Matrix K(6,6);
    double tobs;

    // Measurement loop
    double t = 0;
    double t_old;

    cout << "1";
    for (int i = 1; i <= nobs; i++) {
        //Previous step
        t_old = t;
        Y_old = Y;

        // Time increment and propagation
        Mjd_UTC = obs(i, 1);                       // Modified Julian Date
        t = (Mjd_UTC - Mjd0) * 86400.0;            //Time since epoch [s]

        double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
        IERS(&eopdata, Mjd_UTC, 'l', x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);

        double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
        timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

        double Mjd_TT = Mjd_UTC + TT_UTC / 86400;
        double Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;
        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;

        for (int ii = 1; ii <= 6; ii++) {
            yPhi(ii, 1) = Y_old(ii, 1);
            for (int j = 1; j <= 6; j++) {
                if (ii == j) {
                    yPhi(6 * j + ii, 1) = 1;
                } else {
                    yPhi(6 * j + ii, 1) = 0;
                }
            }
        }

        yPhi = DEInteg(VarEqn, 0, t - t_old, 1e-13, 1e-6, 42, &yPhi);


        // Extract state transition matrices
        for (int j = 1; j <= 6; j++) {
            for (int ii = 1; ii <= 6; ii++) {
                Phi(ii, j) = yPhi(6 * j + ii, 1);
            }
        }

        Y = DEInteg(Accel, 0, t - t_old, 1e-13, 1e-6, 6, &Y_old);

        // Topocentric coordinates
        double theta = gmst(Mjd_UT1);                    // Earth rotation
        Matrix U(3, 3);
        U = R_z(theta);

        Matrix r(3, 1);
        r(1, 1) = Y(1, 1);
        r(2, 1) = Y(2, 1);
        r(3, 1) = Y(3, 1);

        Matrix s(3, 1);
        s = LT * (U * r - Rs);                          // Topocentric position [m]

        // Time update
        TimeUpdate(P, &Phi);

        // Azimuth and partials
        double Azim, Elev;
        Matrix dAds(3, 1);
        Matrix dEds(3, 1);

        AzElPa(&s, Azim, Elev, dAds, dEds);     // Azimuth, Elevation

        Matrix dAdY(6, 1);
        Matrix tempProd(3, 1);
        tempProd = dAds * LT * U;

        dAdY(1, 1) = tempProd(1, 1);
        dAdY(2, 1) = tempProd(2, 1);
        dAdY(3, 1) = tempProd(3, 1);

        dAdY(4, 1) = 0;
        dAdY(5, 1) = 0;
        dAdY(6, 1) = 0;

        // Measurement update
        tobs = obs(i, 2);
        MeasUpdate(tobs, Azim, sigma_az, &dAdY, 6, K, Y, P);


        // Elevation and partials
        r(1, 1) = Y(1, 1);
        r(2, 1) = Y(2, 1);
        r(3, 1) = Y(3, 1);

        s = LT * (U * r - Rs);                          // Topocentric position [m]
        AzElPa(&s, Azim, Elev, dAds, dEds);     // Azimuth, Elevation

        Matrix dEdY(6, 1);
        tempProd = dEds * LT * U;

        dEdY(1, 1) = tempProd(1, 1);
        dEdY(2, 1) = tempProd(2, 1);
        dEdY(3, 1) = tempProd(3, 1);

        dEdY(4, 1) = 0;
        dEdY(5, 1) = 0;
        dEdY(6, 1) = 0;

        // Measurement update
        tobs = obs(i, 3);
        MeasUpdate(tobs, Elev, sigma_el, &dEdY, 6, K, Y, P);

        // Range and partials
        r(1, 1) = Y(1, 1);
        r(2, 1) = Y(2, 1);
        r(3, 1) = Y(3, 1);

        s = LT * (U * r - Rs);                          // Topocentric position [m]

        double Dist = s.norm();

        Matrix dDds(1, 3);

        dDds = (s * (1 / Dist)).transpose();                           // Range


        Matrix dDdY(6, 1);
        tempProd = dDds * LT * U;

        dDdY(1, 1) = tempProd(1, 1);
        dDdY(2, 1) = tempProd(2, 1);
        dDdY(3, 1) = tempProd(3, 1);

        dDdY(4, 1) = 0;
        dDdY(5, 1) = 0;
        dDdY(6, 1) = 0;

        // Measurement update
        tobs = obs(i, 4);
        MeasUpdate(tobs, Dist, sigma_range, &dDdY, 6, K, Y, P);
    }

    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;
    IERS(&eopdata,obs(46,1),'l', x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);

    double UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;
    timediff(UT1_UTC,TAI_UTC, UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);

    double Mjd_TT = Mjd_UTC + TT_UTC/86400;

    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

    Matrix Y0(2,3);

    Y0 = DEInteg (Accel,0,-(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,&Y);

    Matrix Y_true(6,1);
    Y_true(1,1) =5753.173e3;
    Y_true(2,1) = 2673.361e3;
    Y_true(3,1) = 3440.304e3;
    Y_true(4,1) = 4.324207e3;
    Y_true(5,1) = -1.924299e3;
    Y_true(6,1) = -5.728216e3;

    cout << "Error of Position Estimation" << endl;
    cout << "dX " << Y0(1,1)-Y_true(1,1) << "[m]" << endl;
    cout << "dY " << Y0(2,1)-Y_true(2,1) << "[m]" << endl;
    cout << "dZ " << Y0(3,1)-Y_true(3,1) << "[m]" << endl;

    cout << "Error of Velocity Estimation" << endl;
    cout << "dVx " << Y0(4,1)-Y_true(4,1) << "[m]" << endl;
    cout << "dVy " << Y0(5,1)-Y_true(5,1) << "[m]" << endl;
    cout << "dVz " << Y0(6,1)-Y_true(6,1) << "[m]" << endl;

    check_test();
}

