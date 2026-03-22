#include <iostream>
#include <fstream>
#include <cmath>

#include <TFile.h>
#include <TTree.h>

void main(const char* fname){

    //Path must be to a .root file

    TFile *f = TFile::Open(path);

    //TTree initialization

    TTree *t = (TTree*)f->Get("T");

}