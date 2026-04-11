//Standard include

#include <iostream>
#include <fstream>
#include <cmath>
#include <numbers>
#include <chrono>

//ROOT include

#include <TSystem.h> 
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include "RooUnfoldResponse.h"

//Include functions

#include "Functions.h"

//Global variables

//Variable names and axis labels

std::vector<std::string> vars = {"pt", "phi_eta", "y"};

std::vector<std::string> xlabels = {"p^{Z}_{T}[GeV]", "#phi^{*}_{#eta}", "y^{Z}"};

std::vector<std::string> ylabels = {"d#sigma / dp^{Z}_{T}[pb/GeV]", "d#sigma / d#phi^{*}_{#eta} [pb]", "d#sigma / dy^{Z} [pb]"};

//Bounds for drawing histograms;

std::vector<std::pair<Float_t, Float_t>> bounds = {{0, 500}, {0, 15}, {0, 3}};

//Luminosità integrata [pb^{-1}]

double L = 35900;

//Macro to reconstruct Z^0 peak and cross section and obtain histograms of pt, eta and phi;

void CrossSection(const char* fname, std::string outname = "Repo/outFiles/OutDS.root", bool MT = true){

    //Necessary imports for 4-vectors and other utilities, compile macro with + at the end;

    gSystem->Load("libPhysics"); 
    gSystem->Load("libMathCore");
    gSystem->Load("libRooUnfold");
    gROOT->SetBatch(kTRUE);

    //Chrono counter;

    auto start = std::chrono::high_resolution_clock::now();

    //Activating parallel execution by multithreading, if MT is true;

    if (MT){ROOT::EnableImplicitMT();}

    //Initializing DataFrame, fname must be to a .root file;

    ROOT::RDataFrame df("Events", fname);

    /* 

    Z^0 peak is reconstructed with the following filter (based on article):

    1) 2 Muon tracks (nMuon == 2)
    2) Neutral state (Sum of Muon_charge == 0)
    3) Must be isolated (RelIso03 < 0.15)
    4) High transverse momentum (Muon_pt > 25 GeV)
    5) Defined region of Pseudorapidity (|Muon_eta| < 2.4)

    Then a new column is added for invariant mass, selected with abs(Dimuon_mass - 91.1817) < 15;
    We calculate Z^0 transverse momentum, rapidity and optimized angle, we will need them for cross section measurement, definitions are found on article

    */

    auto new_df = df.Filter("nMuon == 2 && "
                              "Muon_charge[0] + Muon_charge[1] == 0 &&"
                              "Muon_pfRelIso03_all[0] < 0.15 && Muon_pfRelIso03_all[1] < 0.15 &&"
                              "Muon_pt[0] > 25 && Muon_pt[1] > 25 &&"
                              "abs(Muon_eta[0]) < 2.4 && abs(Muon_eta[1]) < 2.4")
                      .Define("Z0_p", [](const ROOT::RVecF &pt, const ROOT::RVecF &eta, const ROOT::RVecF &phi, const ROOT::RVecF &mass){
                               ROOT::Math::PtEtaPhiMVector p_1(pt[0],eta[0],phi[0], mass[0]), p_2(pt[1],eta[1], phi[1], mass[1]);
                               return (ROOT::Math::PtEtaPhiMVector)(p_1 + p_2);}, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"})
                      .Define("mass", "Z0_p.M()")
                      .Filter("abs(mass- 91.1817) < 15")
                      .Define("pt", "Z0_p.Pt()")
                      .Define("y", "abs(Z0_p.Rapidity())")
                      .Define("phi_eta", [](const ROOT::RVecF &phi, const ROOT::RVecF &eta){

                                                                float dphi = ROOT::VecOps::DeltaPhi(phi[0], phi[1]);

                                                                float deta = eta[0] - eta[1];

                                                                float sinDteta = 1./std::cosh(deta/2);

                                                                return std::tan((M_PI - dphi)/2.0) * sinDteta;

                                                                }, {"Muon_phi", "Muon_eta"});

    //Output file

    TFile output(outname.c_str(), "RECREATE");

    //Histogram for Invariant Mass;

    auto h_m = new_df.Histo1D({"M_inv", "Z0_mass", 100, 70, 110}, "mass");
    h_m->Draw();
    h_m->Write("InvMass");

    //Loop for diMuon kinematic variables histograms;

    for (int i = 0; i < 3; i++){

        auto h = new_df.Histo1D(ROOT::RDF::TH1DModel(vars[i].c_str(), vars[i].c_str(), 100, bounds[i].first, bounds[i].second), Form("%s", vars[i].c_str()));

        h->Scale(1./L);

        h->Draw("E1 P");
        h->GetXaxis()->SetTitle(xlabels[i].c_str());
        h->GetYaxis()->SetTitle(ylabels[i].c_str());
        h->SetMarkerStyle(20);
        h->SetLineColor(kBlack);

        h->Write(vars[i].c_str());

    }

    //Closing output file

    output.Close();

    //Ending Chrono counter e elapsed time;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    //Elapsed time printing

    std::cout << "Tempo di esecuzione: " << elapsed.count() << " secondi." << std::endl;

    //Deactivating batch mode

    gROOT->SetBatch(kFALSE);

}

//Macro to extract dimuon at generator level, and obtain reco efficiency

void MCSel(const char* fname, std::string outname = "Repo/outFiles/OutMS.root", bool MT = true){

    //Necessary for 4 vectors and other utilities, compile with + at the end

    gSystem->Load("libPhysics");
    gSystem->Load("libMathCore");
    gSystem->Load("libMatrix");
    gSystem->Load("libHist");
    gSystem->Load("libMinuit");
    gSystem->Load("libMathCore");
    gSystem->Load("libRooUnfold");
    gROOT->SetBatch(kTRUE);

    //Chrono counter;

    auto start = std::chrono::high_resolution_clock::now();

    //Optional: activate multithreading

    if (MT){ROOT::EnableImplicitMT();}

    //Initializing DataFrame, fname must be to a .root file;

    ROOT::RDataFrame df("Events", fname);

    //Response matrix

    RooUnfoldResponse response(100, 70, 110, 100, 70, 110);
    
    /*

    Reconstructing Z^0 peak by MC tagging at generator level:

    Loop on GenPart_pdgId, need at least a +13 and -13 con GenPart_pdgId(GenPart_genPartIdxMother) == 23 

    We then compute invariant mass for the pair

    */

    auto new_df = df.Define("genMuon", GenSel, {"nGenPart", "GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass"})
                    .Define("gen_mass", Minv_calculator, {"genMuon"})
                    .Filter(Minv_Range, {"gen_mass"})
                    .Define("gen_pt", Pt_calculator, {"genMuon"})
                    .Define("gen_y", y_calculator, {"genMuon"})
                    .Define("gen_phieta", phi_eta_calculator, {"genMuon"})
                    .Define("recoMuon", Reco, {"nMuon", "Muon_charge", "Muon_pfRelIso03_all","Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"})
                    .Define("reco_mass", Minv_calculator, {"recoMuon"})
                    .Filter(Minv_Range, {"reco_mass"})
                    .Define("reco_pt", Pt_calculator, {"recoMuon"})
                    .Define("reco_y", y_calculator, {"recoMuon"})
                    .Define("reco_phieta", phi_eta_calculator, {"recoMuon"});

    new_df.Foreach([&](Double_t genM, Double_t recoM) {response.Fill((double)recoM, (double)genM);},{"gen_mass", "reco_mass"});

    //Output file

    TFile output(outname.c_str(), "RECREATE");

    //Final histogram
    auto h_mc = new_df.Histo1D({"M_invMC", "DiMuon Mass MC", 100, 70, 110}, "gen_mass");
    h_mc->Draw();
    h_mc->Write("InvMass_MC");

    //Visualization of response matrix

    TH2D* h_response = (TH2D*) response.Hresponse();

    h_response->Draw("COLZ");

    h_response->Write("Response_mass");

    //Closing output file

    output.Close();

    //Ending Chrono counter e elapsed time;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    //Elapsed time printing

    std::cout << "Tempo di esecuzione: " << elapsed.count() << " secondi." << std::endl;

    //Deactivating batch mode

    gROOT->SetBatch(kFALSE);
}