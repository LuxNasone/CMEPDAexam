#include <iostream>
#include <fstream>
#include <cmath>

#include <TFile.h>
#include <TTree.h>

void console(const char* fname){

    //Path must be to a .root file, also a check

    TFile *f = TFile::Open(fname);

    if (!f || f->IsZombie()) { 
        std::cerr << "Error: file not found or corrupted\n"; 
        return;
    }

    //TTree initialization, also a check

    TTree *t = (TTree*)f->Get("Events");

    if (!t) {
        std::cerr << "Error: tree not found in file\n";
        return;
    }

    //Control : plotting Muon_pt histogram

    TH1F *hPT = new TH1F("hPT", "Muon_pt", 100, -50, 200);

    t->Draw("Muon_pt >> hPT");

    hPT->Draw();

}