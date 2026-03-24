#include <iostream>
#include <fstream>
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <Rtypes.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <ROOT/RDataFrame.hxx>

void console(const char* fname){

    //Path must be to a .root file, also a check

    TFile *f = TFile::Open(fname);

    if (!f || f->IsZombie()) { 
        std::cerr << "Error: file not opened\n"; 
        return;
    }

    //TTree initialization, also a check

    TTree *t = (TTree*)f->Get("Events");

    if (!t) {
        std::cerr << "Error: tree not opened\n";
        return;
    }

    //Number of events

    ULong64_t N = t->GetEntries();

    //Muon number, charge, tranverse momentum, pseudorapidity, phi and mass variables

    UInt_t n;
    Int_t q[100];
    Float_t pt[100], eta[100], phi[100], m[100];

    //Address setting

    t->SetBranchAddress("nMuon", &n);
    t->SetBranchAddress("Muon_charge", &q);
    t->SetBranchAddress("Muon_pt", &pt);
    t->SetBranchAddress("Muon_eta", &eta);
    t->SetBranchAddress("Muon_phi", &phi);
    t->SetBranchAddress("Muon_mass", &m);

    //Histogram setup

    TH1D *h = new TH1D("H M_inv", "DiMuon mass", 100, 40, 140);

    //Event loop

    for (ULong64_t k = 0; k < N; k++){

        t->GetEntry(k);

        //Selection : two muon tracks, with opposite charges, high transverse momentum and in certain eta range(article for explanation)

        bool selection = (n == 2 && 
                    q[0] + q[1] == 0 &&
                    pt[0] >= 25 && pt[1] >=25 &&
                    abs(eta[0]) <= 2.4 && abs(eta[1]) <= 2.4);
        
        if (selection){

            //4-impulse of the two muons, with hist fill of invariant mass

            TLorentzVector p_1;
            p_1.SetPtEtaPhiM(pt[0], eta[0], phi[0], m[0]);

            TLorentzVector p_2;
            p_2.SetPtEtaPhiM(pt[1], eta[1], phi[1], m[1]);

            h->Fill((p_1 + p_2).M());

        }
        
    }

    //Resonance plot

    TCanvas *c = new TCanvas("Canvas M_in", "DiMuon canvas", 800, 600);
    h->Draw();
    c->Update();
    c->SaveAs("/home/lux_n/CMEPDA/Exam/Repo/graphs/DiMuonMass.png");

    //Clear

    f->Close();
    delete f;
    delete c;

}