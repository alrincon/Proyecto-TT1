#include "../include/Global.h";
#include "../include/SAT_const.h";
#include <stdio.h>
#include <string.h>

int main(){
    int c = 6;

    extern Matrix eopdata;
    extern Matrix GEOS3data;

    extern Matrix Cnm;
    extern Matrix Snm;

    AUXPARAM AuxParam;

    strcpy(AuxParam.Mjd_UTD, "MJD_UTD");
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;

    FILE *fp1 = fopen("../data/eop19620101.txt", "r");
    if (fp1 == nullptr) {
        printf("error");
        exit(1);
    }

    for (int i = 1; i <= 13; i++) {
        fscanf(fp1, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &((eopdata)(1, i)), &((eopdata)(2, i)), &((eopdata)(3, i)),
               &((eopdata)(4, i)), &((eopdata)(5, i)), &((eopdata)(6, i)),
               &((eopdata)(7, i)), &((eopdata)(8, i)), &((eopdata)(9, i)),
               &((eopdata)(10, i)), &((eopdata)(11, i)), &((eopdata)(12, i)),
               &((eopdata)(13, i)));
    }
    fclose(fp1);

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

        fscanf(fp2, "%d%c%d%c%d %d%c%d%c%lf %lf %lf %lf",
                    year, month, day, hour, min, sec, &((*GEOS3data)(2, i)), &((*GEOS3data)(3, i)), &((*GEOS3data)(4, i))
               );

        &((*GEOS3data)(1, i)) = Mjday(year, month, day, hour,min, sec);
    }
    fclose(fp2);


    FILE *fp3 = fopen("../data/GGM03S.txt", "r");
    if (fp3 == nullptr) {
        printf("error");
        exit(1);
    }

    for (int i = 1; i <= 180; i++) {
        int valor1;
        int valor2;
        double valor5;
        double valor6;



        for(int j = 1; j < i; j++){
            fscanf(fp3, "%d %d %lf %lf %lf %lf",
                                          valor1, valor2, &((*Cnm)(i, j)), &((*Snm)(i, j)), valor5, valor6
            );
        }
    }
    fclose(fp3);



}