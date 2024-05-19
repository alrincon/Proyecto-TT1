#include "../include/Global.h"
#include "../include/SAT_const.h"
#include <stdio.h>
#include <string.h>
#include "../include/Mjday.h"
#include "../include/Position.h"
#include "../include/anglesg.h"


int main(){
    int nobs = 46;
    int infFile = 100;

    Global::initialize();


    /*(*Global::eopdata)(infFile, 13);

    (*Global::obs)(nobs, 4);
    (*Global::Cnm)(181, 181);
    (*Global::Snm)(181, 181);*/

    //fichero GGM03S.txt
    FILE *fp3 = fopen("../data/GGM03S.txt", "r");
    if (fp3 == nullptr) {
        printf("error");
        exit(1);
    }

    for (int i = 0; i < 181; i++) {
        int valor1;
        int valor2;
        double valor5;
        double valor6;



        for(int j = 0; j <= i; j++){
            fscanf(fp3, "%d %d %lf %lf %lf %lf", valor1, valor2, (*Global::Cnm)(i+1, j+1), (*Global::Snm)(i+1, j+1), valor5, valor6);
        }
    }
    fclose(fp3);

    //parametros de orientación de la tierra
    FILE *fp1 = fopen("../data/eop19620101.txt", "r");
    if (fp1 == nullptr) {
        printf("error");
        exit(1);
    }

    for (int i = 1; i <= infFile; i++) {
        fscanf(fp1, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               (*Global::eopdata)(1, i), (*Global::eopdata)(2, i), (*Global::eopdata)(3, i),
               (*Global::eopdata)(4, i), (*Global::eopdata)(5, i), (*Global::eopdata)(6, i),
               (*Global::eopdata)(7, i), (*Global::eopdata)(8, i), (*Global::eopdata)(9, i),
               (*Global::eopdata)(10, i), (*Global::eopdata)(11, i), (*Global::eopdata)(12, i),
               (*Global::eopdata)(13, i));
    }
    fclose(fp1);

    //observaciones
    FILE *fp2 = fopen("../data/GEOS3.txt", "r");
    if (fp2 == nullptr) {
        printf("error");
        exit(1);
    }

    for (int i = 1; i <= 13; i++) {
        int year;
        int month;
        int day;

        int hour;
        int min;
        double sec;

        double az;
        double el;
        double Dist;

        fscanf(fp2, "%d%c%d%c%d %d%c%d%c%lf %lf %lf %lf",
                    year, month, day, hour, min, sec, az, el, Dist);

        (*Global::obs)(1,i) = Mjday(year, month, day, hour,min, sec);
        (*Global::obs)(2,i) = Rad*az;
        (*Global::obs)(3,i) = Rad*el;
        (*Global::obs)(4,i) = 1e3*Dist;
    }
    fclose(fp2);

    double sigma_range = 92.5;          // [m]
    double sigma_az = 0.0224*Rad; // [rad]
    double sigma_el = 0.0139*Rad; // [rad]

    // Kaena Point station
    double lat = Rad*21.5748;     // [rad]
    double lon = Rad*(-158.2706); // [rad]
    double alt = 300.20;                // [m]

    Matrix Rs = Position(lon, lat, alt);

    double Mjd1 = (*Global::obs)(1,1);
    double Mjd2 = (*Global::obs)(9,1);
    double Mjd3 = (*Global::obs)(18,1);

    Matrix r2(1,3);
    Matrix v2(1,3);
    anglesg((*Global::obs)(1,2),(*Global::obs)(9,2),(*Global::obs)(18,2),(*Global::obs)(1,3),(*Global::obs)(9,3),(*Global::obs)(18,3),Mjd1,Mjd2,Mjd3,&Rs,&Rs,&Rs, r2, v2);

    Matrix Y0_apr(2,3);
    Y0_apr(1,1) = r2(1,1); Y0_apr(1,2) = r2(1,2); Y0_apr(1,3) = r2(1,3);
    Y0_apr(2,1) = v2(1,1); Y0_apr(2,2) = v2(1,2); Y0_apr(2,3) = v2(1,3);

    double Mjd0 = Mjday(1995,1,29,02,38,0);

    double Mjd_UTC = (*Global::obs)(9,1);

    Global::AuxParam.Mjd_UTC = Mjd_UTC;
    Global::AuxParam.n = 20;
    Global::AuxParam.m = 20;
    Global::AuxParam.sun = 1;
    Global::AuxParam.moon = 1;
    Global::AuxParam.planets = 1;
    //parametros de modelo

    int n_eqn = 6;


    //añadir mas cosas
}