#include "../include/Global.h"

 void Global::initialize(){
    int nobs = 46;
    int infFile = 100;

    (*Global::eopdata)(infFile, 13);

    (*Global::obs)(nobs, 4);
    (*Global::Cnm)(181, 181);
    (*Global::Snm)(181, 181);
}
