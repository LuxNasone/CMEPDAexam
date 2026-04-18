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

//Number of bins

int n_b = 100;

//Variable names and axis labels

std::vector<std::string> vars = {"pt", "phi_eta", "y"};

std::vector<std::string> xlabels = {"p^{Z}_{T}[GeV]", "#phi^{*}_{#eta}", "y^{Z}"};

std::vector<std::string> ylabels = {"d#sigma / dp^{Z}_{T}[pb/GeV]", "d#sigma / d#phi^{*}_{#eta} [pb]", "d#sigma / dy^{Z} [pb]"};

//Bounds for drawing histograms;

std::vector<std::pair<Float_t, Float_t>> bounds = {{0, 100}, {0, 3}, {0, 1.5}};

//Luminosità integrata [pb^{-1}]

double L = 35900;

//Macro to reconstruct Z^0 peak and cross section and obtain histograms of pt, eta and phi;

std::vector<TH1D> CrossSection(const char* fname,
                                std::string outname = "Repo/outFiles/NotUnfolded.root", 
                                bool MT = true,
                                bool kw = true){

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

    //Loop for histograms;

    std::vector<TH1D> h_v(3);

    for (int i = 0; i < 3; i++){

        auto h_tmp = new_df.Histo1D(ROOT::RDF::TH1DModel(vars[i].c_str(), vars[i].c_str(), n_b, bounds[i].first, bounds[i].second), Form("%s", vars[i].c_str()));

        TH1D h = h_tmp.GetValue(); 

         h_v[i] = h;

        h.Scale(1./L);

    }

    //Saving on file, optional if kw is True

    if(kw){

        TFile output(outname.c_str(), "RECREATE");

        //Histogram for Invariant Mass;

        auto h_m = new_df.Histo1D({"M_inv", "Z0_mass", n_b, 70, 110}, "mass");
        
        h_m->Write("InvMass");

        //Writing for other histograms

        for (size_t i = 0; i < vars.size(); i++){h_v[i].Write(Form("Cross_section(%s)", vars[i].c_str()));}

        //Closing output file

        output.Close();

    }

    //Ending Chrono counter e elapsed time;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    //Elapsed time printing

    std::cout << "Tempo di esecuzione: " << elapsed.count() << " secondi." << std::endl;

    //Deactivating batch mode

    gROOT->SetBatch(kFALSE);

    //return value

    return h_v;

}

std::vector<TH1D> Unfolded(const char* fname,
              const char* rpath = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Response.root", 
              const char* outname = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Unfolded.root"){

    //Matrices file

    TFile* Rfile = TFile::Open(rpath);

    //Response matrix vectors

    std::vector<RooUnfoldResponse> T(vars.size());

    //Filling matrices

    for(size_t i = 0; i < vars.size(); i++){

        TH2D* M = (TH2D*) Rfile->Get(Form("Response_%s", vars[i].c_str()));

        if (!M) {
            std::cerr << "Errore: Response_" << vars[i] << " non trovato!" << std::endl;
            continue;
        }

        TH2D* M_copy = (TH2D*) M->Clone();
                
        RooUnfoldResponse response;

        response.Setup(nullptr, nullptr, M_copy);

        T[i] = response;
    }

    //Closing file

    Rfile->Close();

    //Executing CrossSection

    std::vector<TH1D> h_v = CrossSection(fname, "Repo/outFiles/NotUnfolded.root", true, false);

    TFile* output = new TFile(outname, "RECREATE");

    std::vector<TH1D> h_u(h_v.size());

    //Unfolding

    for (size_t i = 0; i < h_v.size(); i++){

        if (h_v.size() != vars.size()) {std::cerr << "Mismatch size: vars=" << vars.size() << " h_v=" << h_v.size() << std::endl;}

        TH1D* h = (TH1D*)h_v[i].Clone();

        if (!h) {
            std::cerr << "Histogram nullo per i=" << i << std::endl;
            continue;
        }

        RooUnfoldBayes unfold(&T[i], h, 20);

        TH1D* hUnfold = (TH1D*) unfold.Hunfold();

        if (!hUnfold) {
            std::cerr << "Unfold fallito per i=" << i << std::endl;
            continue;
        }

        h_u[i] = *hUnfold;

        hUnfold->Write(Form("%s_unfolded", vars[i].c_str()));

    }

    return h_u;

}

void Comp(const char* fname, 
          const char* outname = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Comparison.root"){

    //Results not-unfolded and unfolded

    std::vector<TH1D> h_nu = CrossSection(fname); 

    std::vector<TH1D> h_u = Unfolded(fname);

    //Outfile

    TFile* output = new TFile(outname, "RECREATE");

    //Plotting 

    if (h_nu.size() == h_u.size()){

        for (size_t i = 0; i < h_u.size(); i++){

            TCanvas* c = new TCanvas(vars[i].c_str(), vars[i].c_str(), 800, 600);

            h_nu[i].SetLineColor(kBlue);
            h_nu[i].SetMarkerColor(kBlue);

            h_nu[i].Draw("E");

            h_u[i].SetLineColor(kRed);
            h_u[i].SetMarkerColor(kRed);

            h_u[i].Draw("E SAME");

            c->Update();

            c->Write();

        }

    }

    output->Close();

}
