#include "../include/Global.h"
#include "../include/SAT_const.h"
#include <stdio.h>
#include <string.h>

/*
int main(){
    extern int nobs;
    extern int infFile;

    extern Matrix eopdata;

    extern Matrix obs;
    extern Matrix Cnm;
    extern Matrix Snm;

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
            fscanf(fp3, "%d %d %lf %lf %lf %lf",
                   valor1, valor2, &(Cnm(i+1, j+1)), &(Snm(i+1, j+1)), valor5, valor6
            );
        }
    }
    fclose(fp3);

    //parametros de modelo
    AUXPARAM AuxParam;
    strcpy(AuxParam.Mjd_UTD, "MJD_UTD");
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;

    //parametros de orientación de la tierra
    FILE *fp1 = fopen("../data/eop19620101.txt", "r");
    if (fp1 == nullptr) {
        printf("error");
        exit(1);
    }

    for (int i = 1; i <= infFile; i++) {
        fscanf(fp1, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &((eopdata)(1, i)), &((eopdata)(2, i)), &((eopdata)(3, i)),
               &((eopdata)(4, i)), &((eopdata)(5, i)), &((eopdata)(6, i)),
               &((eopdata)(7, i)), &((eopdata)(8, i)), &((eopdata)(9, i)),
               &((eopdata)(10, i)), &((eopdata)(11, i)), &((eopdata)(12, i)),
               &((eopdata)(13, i)));
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

        obs.assign(1, i, Mjday(year, month, day, hour,min, sec));
        obs.assign(2,i, Rad*az);
        obs.assign(3,i, Rad*el);
        obs.assign(4,i, 1e3*Dist);
    }
    fclose(fp2);

    double sigma_range = 92.5;          // [m]
    double sigma_az = 0.0224*Rad; // [rad]
    double sigma_el = 0.0139*Rad; // [rad]

    // Kaena Point station
    double lat = Rad*21.5748;     // [rad]
    double lon = Rad*(-158.2706); // [rad]
    double alt = 300.20;                // [m]

    double Rs = Position(lon, lat, alt);

    double Mjd1 = obs(1,1);
    double Mjd2 = obs(9,1);
    double Mjd3 = obs(18,1);

    double r2, v2;
    [r2,v2] = anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);







}*/