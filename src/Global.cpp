#include "../include/Global.h"
#include "../include/Global.h"
#include <fstream>
#include <iostream>
#include <string>



Matrix* Global::eopdata = nullptr;

void Global::eop19620101(int c) {
    Global::eopdata = new Matrix(13, c);
    FILE *fp = fopen("../data/eop19620101.txt", "r");
    if (fp == nullptr) {
        printf("error");
        exit(1);
    }

    for (int i = 1; i <= c; i++) {
        fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &((*Global::eopdata)(1, i)), &((*Global::eopdata)(2, i)), &((*Global::eopdata)(3, i)),
               &((*Global::eopdata)(4, i)), &((*Global::eopdata)(5, i)), &((*Global::eopdata)(6, i)),
               &((*Global::eopdata)(7, i)), &((*Global::eopdata)(8, i)), &((*Global::eopdata)(9, i)),
               &((*Global::eopdata)(10, i)), &((*Global::eopdata)(11, i)), &((*Global::eopdata)(12, i)),
               &((*Global::eopdata)(13, i)));
    }
    fclose(fp);
}

Matrix* Global::GEOS3data = nullptr;

void Global::GEOS3(int c) {
    Global::eopdata = new Matrix(13, c);
    FILE *fp = fopen("../data/GEOS3.txt", "r");
    if (fp == nullptr) {
        printf("error");
        exit(1);
    }

    for (int i = 1; i <= c; i++) {
        string line;

        fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf",
               , ,valor1, valor2, , );
    }
    fclose(fp);
}