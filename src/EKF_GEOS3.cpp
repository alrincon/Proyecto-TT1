#include "../include/Global.h"
#include "../include/SAT_const.h"
#include <stdio.h>
#include <string.h>
#include "../include/Mjday.h"
#include "../include/Position.h"
#include "../include/anglesg.h"


extern Matrix eopdata;
extern Matrix obs;
extern Matrix Cnm;
extern Matrix Snm;
extern aux AuxParam;

int main(){
    int nobs = 46;
    int infFile = 100;

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
            fscanf(fp3, "%d %d %lf %lf %lf %lf", &valor1, &valor2, &Cnm(i+1, j+1), &Snm(i+1, j+1), &valor5, &valor6);
        }
    }

    if (fclose(fp3))   /* Close the stream. */
        perror("fclose error");
    else printf("File mylib/myfile closed successfully.\n");

    //parametros de orientación de la tierra
    FILE *fp1 = fopen("../data/eop19620101.txt", "r");
    if (fp1 == nullptr) {
        printf("error");
        exit(1);
    }

    for (int i = 1; i <= infFile; i++) {
        fscanf(fp1, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &eopdata(i, 1), &eopdata(i, 2), &eopdata(i, 3),
               &eopdata(i, 4), &eopdata(i, 5), &eopdata(i, 6),
               &eopdata(i, 7), &eopdata(i, 8), &eopdata(i, 9),
               &eopdata(i, 10), &eopdata(i, 11), &eopdata(i, 12),
               &eopdata(i, 13));
    }
    if (fclose(fp1))   /* Close the stream. */
        perror("fclose error");
    else printf("File mylib/myfile closed successfully.\n");


    //observaciones
    FILE *fp2 = fopen("../data/GEOS3.txt", "r");
    if (fp2 == nullptr) {
        printf("error");
        exit(1);
    }

    for (int i = 1; i <= nobs; i++) {
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
                    &year, &month, &day, &hour, &min, &sec, &az, &el, &Dist);

        obs(i,1) = Mjday(year, month, day, hour,min, sec);
        obs(i,2) = Rad*az;
        obs(i,3) = Rad*el;
        obs(i,4) = 1e3*Dist;
    }
    if (fclose(fp2))   /* Close the stream. */
        perror("fclose error");
    else printf("File mylib/myfile closed successfully.\n");

    double sigma_range = 92.5;          // [m]
    double sigma_az = 0.0224*Rad; // [rad]
    double sigma_el = 0.0139*Rad; // [rad]

    // Kaena Point station
    double lat = Rad*21.5748;     // [rad]
    double lon = Rad*(-158.2706); // [rad]
    double alt = 300.20;                // [m]

    Matrix Rs = Position(lon, lat, alt);

    double Mjd1 = obs(1,1);
    double Mjd2 = obs(9,1);
    double Mjd3 = obs(18,1);

    cout << "Mjd1 " << Mjd1 << endl;
    cout << "Mjd2 " << Mjd2 << endl;
    cout << "Mjd3 " << Mjd3 << endl;

    Matrix r2(1,3);
    Matrix v2(1,3);

    anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),Mjd1,Mjd2,Mjd3,&Rs,&Rs,&Rs, r2, v2);

    v2.print();

    Matrix Y0_apr(2,3);

    Y0_apr(1,1) = r2(1,1); Y0_apr(1,2) = r2(1,2); Y0_apr(1,3) = r2(1,3);
    Y0_apr(2,1) = v2(1,1); Y0_apr(2,2) = v2(1,2); Y0_apr(2,3) = v2(1,3);

    cout << "funciono" << endl;
    double Mjd0 = Mjday(1995,1,29,02,38,0);

    double Mjd_UTC = obs(9,1);

    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;
    //parametros de modelo

    int n_eqn = 6;


    //añadir mas cosas

    eopdata.print();
}