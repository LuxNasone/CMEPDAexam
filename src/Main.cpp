//Standard include

#include <iostream>
#include <fstream>
#include <cmath>
#include <numbers>
#include <chrono>

//ROOT include

#include <TSystem.h> 
#include <TCanvas.h>
#include <TLegend.h>
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

std::vector<std::pair<Float_t, Float_t>> bounds = {{0, 100}, {0, 3}, {0, 2.5}};

//Range for y axis (purely aesthetic)

std::vector<std::pair<Float_t, Float_t>> range = {{0, 4e4}, {0, 4e4}, {0, 6e3}};

//Luminosità integrata [pb^{-1}]

double L = 35900;

/**
*@brief A raw measurement of distribution for variables used to express differential cross-section (Z trasnverse momentum, rapidity and optimized angle), 
*        without any unfolding procedures appplied.
* This function reads input data from a ROOT file, processes the events, and produces
* three histograms (TH1D) that can be used for further analysis (e.g. unfolding).
*@param fname : file path, required to be a .root file;
*@param outname : name for the output file, required to be a .root file. 
*                  Default : "Repo/outFiles/NotUnfolded.root";
*@param MT :  bool that if true enables multithreading with the option ROOTEnableImplicitMT(). 
*              Default : true;
*@param  mute : bool that if true disables some secondary comments. 
*                Default : false.
*@return a vector of the three histograms (TH1D) saved on the file called outname, meant to be unfolded later on.
*         In order : 
*         [0] : Z tranverse momentum
*         [1] : Optimzed angle = tan(($\pi$ - $\Delta_{\phi}$)/2)/cosh($\Delta_{\eta}$/2)
*         [2] : Z rapidity
*         Advised to visualize with a TBrowser.
*@note Requires ROOT framework and RooUnfold.
*@warning Input file must contain expected tree structure.
*/

std::vector<TH1D> CrossSection(const char* fname,
                                std::string outname = "Repo/outFiles/NotUnfolded.root", 
                                bool MT = true,
                                bool mute = false){

    //Necessary imports for 4-vectors and other utilities, compile macro with + at the end;

    gSystem->Load("libPhysics"); 
    gSystem->Load("libMathCore");
    gSystem->Load("libRooUnfold");
    gROOT->SetBatch(kTRUE);

    //Chrono counter;

    if (mute) {std::cout << "Starting to measure time :" << std::endl;}

    auto start = std::chrono::high_resolution_clock::now();

    //Activating parallel execution by multithreading, if MT is true;

    if (MT){

        ROOT::EnableImplicitMT();

        if (mute) {std::cout << "Activating explicit multithreading, " << "pool size = " << ROOT::GetThreadPoolSize() << std::endl;}

    }

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

    for (size_t i = 0; i < h_v.size(); i++){

        auto h_tmp = new_df.Histo1D(ROOT::RDF::TH1DModel(vars[i].c_str(), vars[i].c_str(), n_b, bounds[i].first, bounds[i].second), Form("%s", vars[i].c_str()));

        TH1D h = h_tmp.GetValue(); 

        h_v[i] = h;

    }

    //Saving on file

    TFile output(outname.c_str(), "RECREATE");

    //Histogram for Invariant Mass;

    auto h_m = new_df.Histo1D({"M_inv", "Z0_mass", n_b, 70, 110}, "mass");
        
    h_m->Write("InvMass");

    //Writing for other histograms

    for (size_t i = 0; i < vars.size(); i++){h_v[i].Write(Form("%s", vars[i].c_str()));}

    //Closing output file

    output.Close();

    //Ending Chrono counter e elapsed time;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    //Elapsed time printing

    if (mute) {std::cout << "Tempo di esecuzione: " << elapsed.count() << " secondi." << std::endl;}

    if (mute) {std::cout << "CrossSection terminato" << std::endl;}

    //Deactivating batch mode

    gROOT->SetBatch(kFALSE);

    //return value

    return h_v;

}

/** 
*@brief Block that applies unfolding on three histograms, with response matrix estimated by Response.cpp (see documentation for more details).
*       Unfolding is applied with RooUnfoldBayes, in the RooUnfold package.
*@param fname : data file path, required to be a .root file;
*@param n_iter : int, number of iteration for unfolding;
*@param rpath : response matrix file path, required to be a .root file;
*@param outname : name for the output file, required to be a .root file. 
*@return A vector of three unfolded histograms, containing same variables analyzed in CrossSection (see documentation for morre details).
*@note Requires ROOT framework and RooUnfold.
*@warning Input file must contain expected tree structure.
*/

std::vector<TH1D> Unfolded(const char* fname,
              int n_iter,
              const char* rpath = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Response.root", 
              const char* outname = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Unfolded.root"){

    //Matrices file

    TFile* Rfile = TFile::Open(rpath);

    if (Rfile->IsZombie()){

        std::cout << "File not opened" << std::endl; 

    }

    //Response matrix vector

    std::vector<RooUnfoldResponse> T(vars.size());

    //Efficiency histogram vector

    std::vector<TH1D> Eff(vars.size());

    //Filling matrices

    for(size_t i = 0; i < vars.size(); i++){

        RooUnfoldResponse* M = (RooUnfoldResponse*) Rfile->Get(Form("Response_%s", vars[i].c_str()));

        TH1D* h_eff = (TH1D*) Rfile->Get(Form("Efficiency for %s", vars[i].c_str()));

        if (!M) {

            std::cerr << "Errore: Response_" << vars[i] << " non trovato!" << std::endl;
            continue;

        }

        T[i] = *M;

        Eff[i] = *h_eff;
    }

    //Closing file

    Rfile->Close();

    //Executing CrossSection

    std::vector<TH1D> h_v = CrossSection(fname, "Repo/outFiles/NotUnfolded.root", true, true);

    std::cout << "This run will be saved on file : " << outname << std::endl;

    TFile* output = new TFile(outname, "RECREATE");

    std::vector<TH1D> h_u(h_v.size());

    //Unfolding

    for (size_t i = 0; i < h_v.size(); i++){

        if (h_v.size() != vars.size()) {
            std::cerr << "Mismatch size: vars =" << vars.size() << " h_v =" << h_v.size() << std::endl;
        }

        TH1D* h = (TH1D*)h_v[i].Clone();

        if (!h) {
            std::cerr << "Istogramma nullo per i=" << i << std::endl;
            continue;
        }

        RooUnfoldBayes unfold(&T[i], h, n_iter);

        TH1D* hUnfold = (TH1D*) unfold.Hunfold();

        if (!hUnfold) {
            std::cerr << "Unfold fallito per i=" << i << std::endl;
            continue;
        }

        /*
        TH1D* hCorr = (TH1D*)hUnfold->Clone("hCorr");

        for (int n = 1; n < n_b; n++){

            double v1 = hUnfold->GetBinContent(n);

            double v2 = Eff[i].GetBinContent(n);

            hCorr->SetBinContent(n, v1 * v2);

        }

        h_u[i] = *hCorr;
        */

        h_u[i] = *hUnfold;

        hUnfold->Write(Form("%s_unfolded", vars[i].c_str()));

    }

    output->Close();

    std::cout << "Saved correctly" << std::endl;

    return h_u;

    std::cout << "Unfolding terminato" << std::endl;

}

/**
*@brief Block that makes comparison between the not unfolded and unfolded distributions listed above.
*@param fname : data file path, required to be a .root file;
*@param n_iter : number of iteration for unfolding;
*@param outname : name for the output file, required to be a .root file. 
*@return nothing, is void type. The graphs can be visualized on a TBrowser
*@note Requires ROOT framework and RooUnfold.
*@warning Input file must contain expected tree structure.
 */

void Comp(const char* fname, 
          int n_iter,
          const char* outname = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Comparison.root"){

    //Results not-unfolded and unfolded

    std::vector<TH1D> h_nu = CrossSection(fname); 

    std::vector<TH1D> h_u = Unfolded(fname, n_iter);

    //Outfile

    if (mute) {std::cout << "This run will be saved on file : " << outname << std::endl;}

    TFile* output = new TFile(outname, "RECREATE");

    if (h_nu.size() == h_u.size()){

        for (size_t i = 0; i < h_u.size(); i++){

            TCanvas* c = new TCanvas(vars[i].c_str(), vars[i].c_str(), 800, 600);

            TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);

            h_nu[i].SetLineColor(kBlue);
            h_nu[i].SetMarkerColor(kBlue);
            h_nu[i].SetMarkerStyle(20);
            h_nu[i].SetStats(0); 

            h_nu[i].GetYaxis()->SetRangeUser(range[i].first, range[i].second);


            legend->AddEntry(&h_nu[i],"Not unfolded", "lep");

            h_nu[i].Draw("E");

            h_u[i].SetLineColor(kRed);
            h_u[i].SetMarkerColor(kRed);
            h_u[i].SetMarkerStyle(20);
            h_u[i].SetStats(0);

            legend->AddEntry(&h_u[i], "Unfolded", "lep");

            h_u[i].Draw("E SAME");

            legend->Draw();

            c->Update();

            c->Write();

            //std::cout << h_nu[i].Integral() << " , " << h_u[i].Integral() << std::endl;

        }

    }

    output->Close();

    if (mute){std::cout << "Saved correctly" << std::endl;}

}
