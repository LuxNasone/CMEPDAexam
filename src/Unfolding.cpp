#include "DataSel.h"
#include "MCSel.h"

//Unfolding procedure, it requires Montecarlo data

void Unfolding(const char* fname){

    //Montecarlo analysis

    MCSel(fname, "Repo/graphs/", true);

    //Data analysis

    DataSel(fname, "Repo/graphs/", true);


}