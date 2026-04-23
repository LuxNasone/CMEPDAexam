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

//Include functions

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "Vars.h"

//Global variables

//Variables name

std::vector<std::string> vars = {"mass", "pt", "phi_eta", "y"};

//Number of bins

int n_b = 100;

//Ranges for variables

std::vector<std::pair<Float_t, Float_t>> bounds = {{0, 100}, {0, 3}, {0, 2.5}};

//Macro to extract dimuon at generator level, and obtain reco efficiency

/**
*@brief Macro to estimate the response matrix, by matching the generated event with reconstructed events and using RooUnfoldResponse
*@param fname : data file path, required to be a .root file;
*@param outname : name for the output file, required to be a .root file.
*                 Default : "Repo/outFiles/Response.root";
*@param MT : bool that if true enables multithreading with the option ROOTEnableImplicitMT(). 
*            Default : true;
*@param mute : bool that if true disables some secondary comments. 
*              Default : false.
*@return nothing, is void type. The graphs can be visualized on a TBrowser
*@note Requires ROOT framework and RooUnfold.
*@warning Input file must contain expected tree structure.
*/

void Response(const char* fname, 
              std::string outname = "Repo/outFiles/Response.root", 
              bool MT = true, 
              bool mute = false){

    //Necessary for 4 vectors and other utilities, compile with + at the end

    gSystem->Load("libPhysics");
    gSystem->Load("libMathCore");
    gSystem->Load("libMatrix");
    gSystem->Load("libHist");
    gSystem->Load("libMinuit");
    gSystem->Load("libMathCore");
    gSystem->AddIncludePath("-I/home/lux_n/CMEPDA/Exam/RooUnfold/src");
    gSystem->Load("/home/lux_n/CMEPDA/Exam/RooUnfold/build/libRooUnfold.so");
    gROOT->SetBatch(kTRUE);

    //Chrono counter;

    if (mute) {std::cout << "Starting to measure time :" << std::endl;}


    auto start = std::chrono::high_resolution_clock::now();

    //Optional: activate multithreading

    if (MT){

        ROOT::EnableImplicitMT();

        if (mute) {std::cout << "Activating explicit multithreading, " << "pool size = " << ROOT::GetThreadPoolSize() << std::endl;}
    
    }

    //Initializing DataFrame, fname must be to a .root file;

    ROOT::RDataFrame df("Events", fname);

    //Response matrix for transverse momentum, opt. angle and rapidity

    std::vector<RooUnfoldResponse> T(bounds.size());

    for (size_t i = 0; i < T.size(); i++){T[i] = RooUnfoldResponse(n_b, bounds[i].first, bounds[i].second, n_b, bounds[i].first, bounds[i].second);}

    /*

    Reconstructing Z^0 peak by MC tagging at generator level:

    Loop on GenPart_pdgId, need at least a +13 and -13 con GenPart_pdgId(GenPart_genPartIdxMother) == 23 

    We then compute invariant mass for the pair

    */

    auto new_df = df.Define("genMuon", GenSel, {"nGenPart", "GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass"})
                    .Define("IsGen", IsTrue, {"nGenPart", "GenPart_pdgId", "GenPart_genPartIdxMother"})
                    .Define("recoMuon", Reco, {"nMuon", "Muon_charge", "Muon_pfRelIso03_all","Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"})
                    .Define("IsReco", IsReco, {"nMuon", "Muon_charge", "Muon_pfRelIso03_all","Muon_pt", "Muon_eta"})
                    .Define("reco_mass", Minv_calculator, {"recoMuon"})
                    .Define("gen_mass", Minv_calculator, {"genMuon"})
                    .Filter(Minv_Range, {"gen_mass"})
                    .Define("gen_pt", Pt_calculator, {"genMuon"})
                    .Define("gen_y", y_calculator, {"genMuon"})
                    .Define("gen_phi_eta", phi_eta_calculator, {"genMuon"})
                    .Define("reco_pt", Pt_calculator, {"recoMuon"})
                    .Define("reco_y", y_calculator, {"recoMuon"})
                    .Define("reco_phi_eta", phi_eta_calculator, {"recoMuon"});

    for (size_t i = 1; i < vars.size(); i++){

        if (vars.size() != (T.size() + 1)){

            std::cout << "Mismatch : T is long " << T.size() << " while vars " << vars.size() << std::endl;

        }

        new_df.Foreach([&](Double_t gen, Double_t reco, bool Gen, bool Reco) {

                            if(Gen && Reco){T[i - 1].Fill((double)reco, (double)gen);}

                            else if (Gen && !Reco){T[i - 1].Miss((double)gen);}

                            else if (!Gen && Reco){T[i - 1].Fake((double)reco);}

                        },{Form("gen_%s", vars[i].c_str()), Form("reco_%s", vars[i].c_str()), "IsGen", "IsReco"});
    }
    
    //Output file

    TFile output(outname.c_str(), "RECREATE");

    //Final histogram
    auto h_mc = new_df.Histo1D({"M_invMC", "DiMuon Mass MC", 100, 70, 110}, "gen_mass");
    h_mc->Write("InvMass_MC");
    auto h_reco = new_df.Histo1D({"M_inv_reco", "DiMuon Mass reco", 100, 70, 110}, "reco_mass");
    h_reco->Write("InvMass_reco");

    //Visualization of response matrix and reconstruction efficiency

    for (size_t i = 1; i < vars.size(); i++){

        T[i - 1].Write(Form("Response_%s", vars[i].c_str()));

        TH2D* M = (TH2D*) T[i - 1].Hresponse();

        TH1D* h_true = (TH1D*) T[i - 1].Htruth();

        TH1D* h_matched = M->ProjectionY();

        TH1D* h_eff = (TH1D*) h_matched->Clone();

        h_eff->Divide(h_matched, h_true, 1.0, 1.0, "B");

        h_eff->Write(Form("Efficiency for %s", vars[i].c_str()));

    }

    //Validation for response matrix

    for (size_t i = 1; i < vars.size(); i++){

        if (vars.size() != (bounds.size() + 1)){

            std::cout << "Mismatch : bounds is long " << bounds.size() << " while vars " << vars.size() << std::endl;

        }


        TCanvas* c = new TCanvas(vars[i].c_str(), vars[i].c_str(), 800, 600);

        auto h_true_ptr = new_df.Histo1D({Form("%s_MC", vars[i].c_str()), Form("%s_MC", vars[i].c_str()), 100, bounds[i - 1].first, bounds[i - 1].second}, Form("gen_%s", vars[i].c_str()));

        auto h_obt_ptr = new_df.Histo1D({Form("%s_Reco", vars[i].c_str()), Form("%s_Reco", vars[i].c_str()), 100, bounds[i - 1].first, bounds[i - 1].second}, Form("reco_%s", vars[i].c_str()));

        TH1D h_true = *h_true_ptr;
        TH1D h_obt = *h_obt_ptr;

        RooUnfoldBayes unfold(&T[i-1], &h_obt, 20);

        TH1D* hUnfold = (TH1D*) unfold.Hunfold();

        hUnfold->SetLineColor(kBlue);
        hUnfold->SetMarkerColor(kBlue);
        hUnfold->SetMarkerStyle(20);

        hUnfold->Draw("P E");

        h_true.Draw("SAME");

        c->Update();

        c->Write(Form("%s_unfolded", vars[i].c_str()));


    }

    //Closing output file

    output.Close();

    //Ending Chrono counter e elapsed time;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    //Elapsed time printing

    if (mute) {std::cout << "Tempo di esecuzione: " << elapsed.count() << " secondi." << std::endl;}

    if (mute) {std::cout << "Response terminato" << std::endl;}

    //Deactivating batch mode

    gROOT->SetBatch(kFALSE);
}