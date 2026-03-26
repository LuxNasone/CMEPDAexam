#include <iostream>
#include <fstream>
#include <cmath>

#include <TSystem.h> 
#include <TFile.h>
#include <TTree.h>
#include <Rtypes.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <ROOT/RDataFrame.hxx>

void Main(const char* fname){

    //Necessary imports for 4-vectors, compile macro with + at the end

    gSystem->Load("libPhysics"); 

    //Import necessary to save graphs

    std::string outdir = std::string(getenv("HOME")) + "/CMEPDA/Exam/Repo/graphs/";

    gSystem->mkdir(outdir.c_str(), true);

    //Initializing DataFrame, fname must be to a .root file

    ROOT::RDataFrame df("Events", fname);

    /* 

    Z^0 peak is reconstructed with the following filter (based on article):

    1) 2 Muon tracks (nMuon == 2)
    2) Neutral state (Sum of Muon_charge == 0)
    3) Must be isolated (RelIso03 < 0.15)
    4) High transverse momentum (Muon_pt > 25 GeV)
    5) Defined region of Pseudorapidity (|Muon_eta| < 2.4)

    */

    auto filtered = df.Filter("nMuon == 2 && "
                              "Muon_charge[0] + Muon_charge[1] == 0 &&"
                              "Muon_pfRelIso03_all[0] < 0.15 && Muon_pfRelIso03_all[1] < 0.15 &&"
                              "Muon_pt[0] > 25 && Muon_pt[1] > 25 &&"
                              "abs(Muon_eta[0]) < 2.4 && abs(Muon_eta[1]) < 2.4");

    //Need to define a new column with invariant mass, also restrict to value near Z^0 mass : abs(M_in - 91.1817) < 15

    auto df2 = filtered.Define("Dimuon_mass", [](const ROOT::RVecF &pt, const ROOT::RVecF &eta, const ROOT::RVecF &phi, const ROOT::RVecF &mass){
                               ROOT::Math::PtEtaPhiMVector p_1(pt[0],eta[0],phi[0], mass[0]);
                               ROOT::Math::PtEtaPhiMVector p_2(pt[1],eta[1],phi[1], mass[1]);
                               return (Float_t)(p_1 + p_2).M();}, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"});

    auto df3 = df2.Filter("abs(Dimuon_mass - 91.1817) < 15");

    //Histogram for Invariant Mass

    auto h_dmm = df3.Histo1D({"M_inv", "DiMuon_mass", 100, 70, 110}, "Dimuon_mass");
    TCanvas *c_dmm = new TCanvas("M_inv", "DiMuon_mass_canvas", 800, 600);
    h_dmm->Draw();
    c_dmm->Update();
    c_dmm->SaveAs((outdir + "DiMuonMass.png").c_str());

    //Histogram for transverse momentum

    auto h_pt = df3.Histo1D({"Pt", "Transverse_momentum", 100, 25, 200}, "Muon_pt");
    TCanvas *c_pt = new TCanvas("Pt", "Transverse_momentum_canvas", 800, 600);
    h_pt->Draw();
    c_pt->Update();
    c_pt->SaveAs((outdir +"DiMuonPT.png").c_str());

    //Histogram for pseudorapidity

    auto h_eta = df3.Histo1D({"Eta", "Muon_pseudorapidity", 100, -2.4, 2.4}, "Muon_eta");
    TCanvas *c_eta = new TCanvas("Eta", "Pseudorapidity_canvas", 800, 600);
    h_eta->Draw();
    c_eta->Update();
    c_eta->SaveAs((outdir + "DiMuonEta.png").c_str());

    //Histogram for azimuthal angle

    auto h_phi = df3.Histo1D({"Phi", "Azimuthal_angle", 100, -3.14, 3.14}, "Muon_phi");
    TCanvas *c_phi = new TCanvas("Phi", "Azimuthal_angle_canvas", 800, 600);
    h_phi->Draw();
    c_phi->Update();
    c_phi->SaveAs((outdir + "DiMuonPhi.png").c_str());

    /*
    //Histogram setup

    TH1D *h = new TH1D("H M_inv", "DiMuon mass", 100, 70, 110);

    //Event (explicit) loop

    for (ULong64_t k = 0; k < N; k++){

        t->GetEntry(k);

        //Selection : two muon tracks, with opposite charges, high transverse momentum and in certain eta range(article for explanation)

        bool selection = ((n == 2) && (q[0] + q[1] == 0) && (pt[0] >= 25 && pt[1] >=25) && (abs(eta[0]) <= 2.4 && abs(eta[1]) <= 2.4) && (iso[0] < 0.15 && iso[1] < 0.15));
        
        if (selection) {

            //4-impulse of the two muons, with hist fill of invariant mass

            TLorentzVector p_1, p_2;

            p_1.SetPtEtaPhiM(pt[0], eta[0], phi[0], m[0]);
            p_2.SetPtEtaPhiM(pt[1], eta[1], phi[1], m[1]);

            Float_t MdM = (p_1 + p_2).M();

            if (abs(MdM - 91.1876) < 15){ h->Fill(MdM); }

        }
        
    }

    //Resonance plot

    TCanvas *c = new TCanvas("Canvas M_in", "DiMuon canvas", 800, 600);
    h->Draw();
    c->Update();
    c->SaveAs("/home/lux_n/CMEPDA/Exam/Repo/graphs/DiMuonMass.png");

    //Clear

    */

    //delete c;

}