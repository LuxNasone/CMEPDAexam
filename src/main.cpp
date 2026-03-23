#include <iostream>
#include <fstream>
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <Rtypes.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLeaf.h>
#include <ROOT/RDataFrame.hxx>

Float_t M_inv(Float_t pt[100], Float_t eta[100], Float_t m[100]){

    Float_t pT = pt[0] + pt[1];

    Float_t pl = pt[0] * std::sinh(eta[0]) + pt[1] * std::sinh(eta[1]);

    Float_t p_sq = pT*pT + pl*pl;

    Float_t E = std::sqrt(pt[0] * pt[0] * cosh(eta[0]) * cosh(eta[0]) + m[0] * m[0]) + std::sqrt(pt[1] * pt[1] * cosh(eta[1]) * cosh(eta[1]) + m[1] * m[1]); 

    return std::sqrt(E*E - p_sq);

}

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

    //Number of events

    ULong64_t N = t->GetEntries();

    //Branch for muon number

    UInt_t n;
    t->SetBranchAddress("nMuon", &n);

    //Branch for muon charges

    Int_t q[100];
    t->SetBranchAddress("Muon_charge", &q);

    //Branch for transverse momentum

    Float_t pt[100];
    t->SetBranchAddress("Muon_pt", &pt);

    //Branch for pseudorapdity

    Float_t eta[100];
    t->SetBranchAddress("Muon_eta", &eta);

    //Branch for muon masses

    Float_t m[100];
    t->SetBranchAddress("Muon_mass", &m);

    //Histogram setup

    TH1D *h = new TH1D("H M_inv", "DiMuon mass", 100, 0, 100);

    //Event selection

    for (ULong64_t k = 0; k < N; k++){

        t->GetEntry(k);

        
        bool selection = (n == 2 && 
                    q[0] + q[1] == 0 &&
                    pt[0] >= 25 && pt[1] >=25 &&
                    abs(eta[0]) <= 2.4 && abs(eta[1]) <= 2.4);
        
        if (selection){

            h->Fill(M_inv(pt, eta, m));

        }
        
    }

    TCanvas *c = new TCanvas("Canvas M_in", "DiMuon canvas", 800, 600);
    h->Draw();
    c->Update();
    c->SaveAs("DiMuonMass.png");

    f->Close();
    delete f;
    delete c;

}