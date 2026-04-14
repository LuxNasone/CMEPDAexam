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
#include "RooUnfoldBayes.h"

//Global variables

//Variable names and axis labels

std::vector<std::string> vars = {"pt", "phi_eta", "y"};

std::vector<std::string> xlabels = {"p^{Z}_{T}[GeV]", "#phi^{*}_{#eta}", "y^{Z}"};

std::vector<std::string> ylabels = {"d#sigma / dp^{Z}_{T}[pb/GeV]", "d#sigma / d#phi^{*}_{#eta} [pb]", "d#sigma / dy^{Z} [pb]"};

//Bounds for drawing histograms;

std::vector<std::pair<Float_t, Float_t>> bounds = {{0, 100}, {0, 15}, {0, 2.5}};

//Luminosità integrata [pb^{-1}]

double L = 35900;

//Macro to reconstruct Z^0 peak and cross section and obtain histograms of pt, eta and phi;

void CrossSection(const char* fname, 
                  const char* rpath = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Response.root" ,
                  std::string outname = "Repo/outFiles/Out.root", 
                  bool MT = true){

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

    //Opening file with response matrices

    TFile* Rfile = TFile::Open(rpath);

    //Getting TH2Ds

    TH2D* Response_pt  = (TH2D*)Rfile->Get("Response_pt");
    TH2D* Response_phieta  = (TH2D*)Rfile->Get("Response_phieta");
    TH2D* Response_y  = (TH2D*)Rfile->Get("Response_y");

    //Initializing Response matrices

    RooUnfoldResponse T_pt(100, 0, 100, 100, 0, 100);
    RooUnfoldResponse T_phieta(100, 0, 15, 100, 0, 15);
    RooUnfoldResponse T_y(100, 0, 2.5, 100, 0, 2.5);

    for (int i = 1; i <= 100; i++){

        for (int j = 1; j <= 100; j++){

            double weight_pt = Response_pt->GetBinContent(i,j);

            T_pt.Fill(j-1, i-1, weight_pt);

            double weight_phieta = Response_phieta->GetBinContent(i,j);

            T_phieta.Fill(j-1, i-1, weight_phieta);

            double weight_y = Response_y->GetBinContent(i,j);

            T_y.Fill(j-1, i-1, weight_y);

        }
    }

    std::vector<RooUnfoldResponse> T = {T_pt, T_phieta, T_y};

    //Closing file

    Rfile->Close();

    //Output file

    TFile output(outname.c_str(), "RECREATE");

    //Histogram for Invariant Mass;

    auto h_m = new_df.Histo1D({"M_inv", "Z0_mass", 100, 70, 110}, "mass");
    h_m->Draw();
    h_m->Write("InvMass");

    //Loop for diMuon kinematic variables histograms;

    for (int i = 0; i < 3; i++){

        auto h_tmp = new_df.Histo1D(ROOT::RDF::TH1DModel(vars[i].c_str(), vars[i].c_str(), 100, bounds[i].first, bounds[i].second), Form("%s", vars[i].c_str()));

        TH1D* h = h_tmp.GetPtr();
        RooUnfoldBayes unfold(&T[i], h, 4);

        TH1D* hUnfold = (TH1D*) unfold.Hunfold();

        hUnfold->Draw("E1 P");
        hUnfold->GetXaxis()->SetTitle(xlabels[i].c_str());
        hUnfold->GetYaxis()->SetTitle(ylabels[i].c_str());
        hUnfold->SetMarkerStyle(20);
        hUnfold->SetLineColor(kBlack);

        hUnfold->Write(vars[i].c_str());

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
