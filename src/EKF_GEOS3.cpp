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
#include "../include/largestRoot.h"


#include "../include/EFK_Tests.h"


extern Matrix eopdata;
extern Matrix obs;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix PC;
extern aux AuxParam;


int main(){
    int nobs = 46;
    int infFile = 21413;

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
               &eopdata( 1,i), &eopdata( 2,i), &eopdata( 3,i),
               &eopdata(4,i), &eopdata(5,i), &eopdata(6,i),
               &eopdata(7,i), &eopdata(8,i), &eopdata(9,i),
               &eopdata(10,i), &eopdata(11,i), &eopdata(12,i),
               &eopdata(13,i));
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

    Matrix r2(3,1);
    Matrix v2(3,1);

    anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),Mjd1,Mjd2,Mjd3,&Rs,&Rs,&Rs, r2, v2);
    Matrix Y0_apr(1,6);

    Y0_apr(1,1) = r2(1,1);
    Y0_apr(1,2) = r2(2,1);
    Y0_apr(1,3) = r2(3,1);
    Y0_apr(1,4) = v2(1,1);
    Y0_apr(1,5) = v2(2,1);
    Y0_apr(1,6) = v2(3,1);


    double Mjd0 = Mjday(1995,1,29,02,38,0);

    double Mjd_UTC = obs(9,1);

    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_UTC;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;
    //parametros de modelo

    int n_eqn  = 6;

    Matrix Y(1,6);
    Matrix Y_old(1,6);

    Y0_apr.print();
    cout << "Resultado" << endl;
    Y = DEInteg(Accel,0.0,-(obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,&Y0_apr);
    Y.print();

    Matrix P(6,6);

    for (int i = 1; i <= 3; i++) {
        P(i, i) = 1e8;
    }

    for (int i = 4; i <= 6; i++) {
        P(i, i) = 1e3;
    }

    Matrix LT(3,3);
    LT = LTC(lon,lat);

    Matrix yPhi(1,42);
    Matrix Phi(6,6);

    Matrix K(6,6);
    double tobs;

    // Measurement loop
    double t = 0;
    double t_old;

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
            yPhi(1,ii) = Y_old(1,ii);
            for (int j = 1; j <= 6; j++) {
                if (ii == j) {
                    yPhi(1,6 * j + ii) = 1;
                } else {
                    yPhi(1,6 * j + ii) = 0;
                }
            }
        }

        yPhi = DEInteg(VarEqn, 0.0, t - t_old, 1e-13, 1e-6, 42, &yPhi);

        // Extract state transition matrices
        for (int j = 1; j <= 6; j++) {
            for (int ii = 1; ii <= 6; ii++) {
                Phi(j,ii) = yPhi(1,6 * j + ii);
            }
        }

        Y = DEInteg(Accel, 0.0, t - t_old, 1e-13, 1e-6, 6, &Y_old);
        // Topocentric coordinates
        double theta = gmst(Mjd_UT1);                    // Earth rotation
        Matrix U(3, 3);
        U = R_z(theta);

        Matrix r(3, 1);
        r(1, 1) = Y(1, 1);
        r(2, 1) = Y(1, 2);
        r(3, 1) = Y(1, 3);

        Matrix s(3, 1);
        s = LT * (U * r - Rs);                          // Topocentric position [m]

        // Time update
        TimeUpdate(P, &Phi);

        // Azimuth and partials
        double Azim, Elev;
        Matrix dAds(3, 1);
        Matrix dEds(3, 1);

        AzElPa(&s, Azim, Elev, dAds, dEds);     // Azimuth, Elevation

        Matrix dAdY(1, 6);
        Matrix tempProd(1, 3);
        tempProd = dAds.transpose() * LT * U;

        dAdY(1, 1) = tempProd(1, 1);
        dAdY(1, 2) = tempProd(1, 2);
        dAdY(1, 3) = tempProd(1, 3);

        dAdY(1, 4) = 0;
        dAdY(1, 5) = 0;
        dAdY(1, 6) = 0;

        // Measurement update
        tobs = obs(i, 2);

        //falla aqui
        MeasUpdate(tobs, Azim, sigma_az, &dAdY, 6, K, Y, P);


        // Elevation and partials
        r(1, 1) = Y(1, 1);
        r(2, 1) = Y(1, 2);
        r(3, 1) = Y(1, 3);

        s = LT * (U * r - Rs);                          // Topocentric position [m]
        AzElPa(&s, Azim, Elev, dAds, dEds);     // Azimuth, Elevation

        Matrix dEdY(1, 6);
        tempProd = dEds.transpose() * LT * U;

        dEdY(1, 1) = tempProd(1, 1);
        dEdY(1, 2) = tempProd(1, 2);
        dEdY(1, 3) = tempProd(1, 3);

        dEdY(1, 4) = 0;
        dEdY(1, 5) = 0;
        dEdY(1, 6) = 0;

        // Measurement update
        tobs = obs(i, 3);
        MeasUpdate(tobs, Elev, sigma_el, &dEdY, 6, K, Y, P);

        // Range and partials
        r(1, 1) = Y(1, 1);
        r(2, 1) = Y(1, 2);
        r(3, 1) = Y(1, 3);

        s = LT * (U * r - Rs);                          // Topocentric position [m]

        double Dist = s.norm();

        Matrix dDds(1, 3);

        dDds = (s * (1 / Dist)).transpose();                           // Range


        Matrix dDdY(1, 6);
        tempProd = dDds * LT * U;

        dDdY(1, 1) = tempProd(1, 1);
        dDdY(1, 2) = tempProd(1, 2);
        dDdY(1, 3) = tempProd(1, 3);

        dDdY(1, 4) = 0;
        dDdY(1, 5) = 0;
        dDdY(1, 6) = 0;

        // Measurement update
        tobs = obs(i, 4);
        MeasUpdate(tobs, Dist, sigma_range, &dDdY, 6, K, Y, P);

        cout << "Ciclo " << i << endl;
    }
    cout << "Fin ciclos" << endl;

    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;
    IERS(&eopdata,obs(46,1),'l', x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);

    double UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;
    timediff(UT1_UTC,TAI_UTC, UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);

    double Mjd_TT = Mjd_UTC + TT_UTC/86400;

    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

    Matrix Y0(1,6);

    Y0 = DEInteg (Accel,0,-(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,&Y);

    Matrix Y_true(1,6);
    Y_true(1,1) =5753.173e3;
    Y_true(1,2) = 2673.361e3;
    Y_true(1,3) = 3440.304e3;
    Y_true(1,4) = 4.324207e3;
    Y_true(1,5) = -1.924299e3;
    Y_true(1,6) = -5.728216e3;

    cout << "Error of Position Estimation" << endl;
    cout << "dX " << Y0(1,1)-Y_true(1,1) << "[m]" << endl;
    cout << "dY " << Y0(1,2)-Y_true(1,2) << "[m]" << endl;
    cout << "dZ " << Y0(1,3)-Y_true(1,3) << "[m]" << endl;

    cout << "Error of Velocity Estimation" << endl;
    cout << "dVx " << Y0(1,4)-Y_true(1,4) << "[m]" << endl;
    cout << "dVy " << Y0(1,5)-Y_true(1,5) << "[m]" << endl;
    cout << "dVz " << Y0(1,6)-Y_true(1,6) << "[m]" << endl;

    //check_test();
}

