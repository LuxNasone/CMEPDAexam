#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

#include <TSystem.h> 
#include <TFile.h>
#include <TTree.h>
#include <Rtypes.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <ROOT/RDataFrame.hxx>

void Main(const char* fname, std::string outname){

    //Necessary imports for 4-vectors, compile macro with + at the end

    gSystem->Load("libPhysics"); 

    //Chrono counter

    auto start = std::chrono::high_resolution_clock::now();

    //Activating parallel execution by multithreading (kind of a black box ... maybe put it as an option?)
    ROOT::EnableImplicitMT();

    //Initializing DataFrame, fname must be to a .root file

    ROOT::RDataFrame df("Events", fname);

    /* 

    Z^0 peak is reconstructed with the following filter (based on article):

    1) 2 Muon tracks (nMuon == 2)
    2) Neutral state (Sum of Muon_charge == 0)
    3) Must be isolated (RelIso03 < 0.15)
    4) High transverse momentum (Muon_pt > 25 GeV)
    5) Defined region of Pseudorapidity (|Muon_eta| < 2.4)

    Then a new column is added for invariant mass, selected with abs(Dimuon_mass - 91.1817) < 15

    */

    auto new_df = df.Filter("nMuon == 2 && "
                              "Muon_charge[0] + Muon_charge[1] == 0 &&"
                              "Muon_pfRelIso03_all[0] < 0.15 && Muon_pfRelIso03_all[1] < 0.15 &&"
                              "Muon_pt[0] > 25 && Muon_pt[1] > 25 &&"
                              "abs(Muon_eta[0]) < 2.4 && abs(Muon_eta[1]) < 2.4")
                      .Define("Dimuon_mass", [](const ROOT::RVecF &pt, const ROOT::RVecF &eta, const ROOT::RVecF &phi, const ROOT::RVecF &mass){
                               ROOT::Math::PtEtaPhiMVector p_1(pt[0],eta[0],phi[0], mass[0]), p_2(pt[1],eta[1], phi[1], mass[1]);
                               return (Float_t)(p_1 + p_2).M();}, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"})
                      .Filter("abs(Dimuon_mass - 91.1817) < 15")
                      .Define("Muon0_pt",  "Muon_pt[0]")
                      .Define("Muon1_pt",  "Muon_pt[1]")
                      .Define("Muon0_eta", "Muon_eta[0]")
                      .Define("Muon1_eta", "Muon_eta[1]")
                      .Define("Muon0_phi", "Muon_phi[0]")
                      .Define("Muon1_phi", "Muon_phi[1]");

    //Histogram for Invariant Mass

    auto h_dmm = new_df.Histo1D({"M_inv", "DiMuon_mass", 100, 70, 110}, "Dimuon_mass");
    TCanvas *c_dmm = new TCanvas("M_inv", "DiMuon_mass_canvas", 800, 600);
    h_dmm->Draw();
    c_dmm->Update();
    c_dmm->SaveAs((outname+"Mass.png").c_str());

    delete c_dmm;

    //Variables necessary for loop

    std::vector<std::string> vars = {"pt", "eta", "phi"};
    std::vector<std::pair<Float_t, Float_t>> bounds = {{15, 200},
                                              {-3, 3},
                                              {-4, 4}};

    //Loop for diMuon kinematic variables histogram

    for (int i = 0; i < 3; i++){

        TCanvas *c = new TCanvas((vars[i] + "canvas").c_str(), vars[i].c_str(), 800, 600);

        auto h_1 = new_df.Histo1D(ROOT::RDF::TH1DModel(vars[i].c_str(), (vars[i] + to_string(1)).c_str(), 100 ,bounds[i].first, bounds[i].second), 
                                  ("Muon" + to_string(0) + "_" + vars[i]).c_str());
        auto h_2 = new_df.Histo1D(ROOT::RDF::TH1DModel(vars[i].c_str(), (vars[i] + to_string(2)).c_str(), 100 , bounds[i].first, bounds[i].second), 
                                  ("Muon" + to_string(1) + "_" + vars[i]).c_str());

        h_1->SetLineColor(kBlue);
        h_2->SetLineColor(kRed);

        h_1->Draw();
        h_2->Draw("SAME");

        c->Update();
        c->SaveAs((outname + vars[i] + ".png").c_str());

        delete c;
    }


    //Ending Chrono counter e elapsed time

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    //Print del tempo

    std::cout << "Tempo di esecuzione: " << elapsed.count() << " secondi." << std::endl;

}