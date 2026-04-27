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

//Unfolding toolbox

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

//Useful functions

#include "Vars.h"

/** @defgroup GlobalVariables Global Variables */

/// @ingroup GlobalVariables
/// @brief Number of bins, used both for histograms and response matrix 

int n_b = 100;

/// @ingroup GlobalVariables
/// @brief variable names, used in loops over RDataFrame columns 

std::vector<const char*> vars = {"pt", "phi_eta", "y"};

/// @ingroup GlobalVariables
/// @brief x axis names, used in loops for plots. These are type std::string, must be converted to char to be used in some application, using c_str() method 

std::vector<const char*> xlabels = {"p^{Z}_{T}[GeV]", "#phi^{*}_{#eta}", "y^{Z}"};

/// @ingroup GlobalVariables
/// @brief x axis names, used in loops for plots. These are type std::string, must be converted to char to be used in some application, using c_str() method


std::vector<const char*> ylabels = {"d#sigma / dp^{Z}_{T}[pb/GeV]", "d#sigma / d#phi^{*}_{#eta} [pb]", "d#sigma / dy^{Z} [pb]"};

/// @ingroup GlobalVariables
/// @brief bounds for variables, both for graphs but also for ranges in response matrix estimation 

std::vector<std::pair<Float_t, Float_t>> bounds = {{0, 100}, {0, 3}, {0, 2.5}};

/// @ingroup GlobalVariables
/// @brief bounds for y axis, purely aesthetic 

std::vector<std::pair<Float_t, Float_t>> range = {{0, 8.5e5}, {0, 1.2e6}, {0, 1.2e5}};

/// @ingroup GlobalVariables
/// @brief integrated luminosity for data used, expressed in [pb^{-1}]

double L = 8740.119304;

/// @ingroup GlobalVariables
/// @brief sigma for Z production in two muons [pb]

double s = 204;

/**
*@brief Macro to reconstruct the distributions of variables used to express the differential cross-section (Z transverse momentum, rapidity and optimized angle), 
*        without any unfolding procedures applied.  
* This function reads input data from a ROOT file in NANOAD format, processes the events, and produces
* three histograms (TH1D) that can be used for further analysis (e.g. unfolding).
*@param folder_name : folder path, required to contain at least a proper .root file;
*@param outname : name for the output file, advised to use a .root file in order to easily see results with a TBrowser. 
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
*@note Requires ROOT framework. 
*@warning Input file must have expected structure, like NanoAD CMS OpenData.
*/

std::vector<TH1D> CrossSection(const char* folder_name,
                                const char* outname = "Repo/outFiles/NotUnfolded.root", 
                                bool MT = true,
                                bool mute = false){

    //Necessary imports for 4-vectors and other utilities, compile macro with + at the end;

    gSystem->Load("libPhysics"); 
    gSystem->Load("libMathCore");
    gSystem->Load("libRooUnfold");
    gROOT->SetBatch(kTRUE);

    //Chrono counter;

    if (!mute) {std::cout << "Starting to measure time :" << std::endl;}

    auto start = std::chrono::high_resolution_clock::now();

    //Activating parallel execution by multithreading, if MT is true;

    if (MT){

        ROOT::EnableImplicitMT();

        if (!mute) {std::cout << "Activating explicit multithreading, " << "pool size = " << ROOT::GetThreadPoolSize() << std::endl;}

    }

    //Initializing DataFrame, folder_name must contain a .root file;

    ROOT::RDataFrame df("Events", Form("%s/*",folder_name));

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

        auto h_tmp = new_df.Histo1D(ROOT::RDF::TH1DModel(vars[i], vars[i], n_b, bounds[i].first, bounds[i].second), Form("%s", vars[i]));

        TH1D h = h_tmp.GetValue(); 

        h_v[i] = h;

    }

    //Saving on file

    TFile output(outname, "RECREATE");

    //Histogram for Invariant Mass;

    auto h_m = new_df.Histo1D({"M_inv", "Z0_mass", n_b, 70, 110}, "mass");
        
    h_m->Write("InvMass");

    //Writing for other histograms

    for (size_t i = 0; i < vars.size(); i++){h_v[i].Write(Form("%s", vars[i]));}

    //Closing output file

    output.Close();

    //Ending Chrono counter e elapsed time;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    //Elapsed time printing

    if (!mute) {std::cout << "Tempo di esecuzione: " << elapsed.count() << " secondi." << std::endl;}

    if (!mute) {std::cout << "CrossSection terminato" << std::endl;}

    //Deactivating batch mode

    gROOT->SetBatch(kFALSE);

    //return value

    return h_v;

}

//Macro to extract dimuon at generator level, and obtain reco efficiency

/**
 * @brief Calculates the  response matrices by matching generated and reconstructed events of interest (Z in dimuon).
 *
 *        The function builds three RooUnfoldResponse objects (one per observable)
 *        and stores them in an output ROOT file.
 *
 * @param folder_name Path to the input folder cointaining ROOT file with generated and reconstructed data.
 *              The file must follow the expected tree structure of the NANOAD format.
 * @param outname Path to the output ROOT file where the response matrices are saved.
 *                Default: "Repo/outFiles/Response.root".
 * @param MT If true, enables ROOT implicit multithreading via ROOT::EnableImplicitMT().
 *           This improves performance but does not affect results.
 *           Default: true.
 * @param mute If true, suppresses non-essential log messages printed to stdout.
 *             Default: false.
 *
 * @return void. The function writes RooUnfoldResponse objects and related histograms
 *         to the output ROOT file, accessible via TBrowser.
 *
 * @note Requires ROOT framework and RooUnfold package.
 *
 * @warning Input file must have expected structure, like NanoAD CMS OpenData.
 */


void Response(const char* folder_name, 
              const char* outname = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Response.root", 
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

    if (!mute) {std::cout << "Starting to measure time :" << std::endl;}


    auto start = std::chrono::high_resolution_clock::now();

    //Optional: activate multithreading

    if (MT){

        ROOT::EnableImplicitMT();

        if (!mute) {std::cout << "Activating explicit multithreading, " << "pool size = " << ROOT::GetThreadPoolSize() << std::endl;}
    
    }

    //Initializing DataFrame, folder must contain a .root file;

    ROOT::RDataFrame df("Events", Form("%s",folder_name));

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

    for (size_t i = 0; i < vars.size(); i++){

        if (vars.size() != (T.size())){

            std::cerr << "Mismatch : T is long " << T.size() << " while vars " << vars.size() << std::endl;

        }

        new_df.Foreach([&](Double_t gen, Double_t reco, bool Gen, bool Reco, Float_t w_1, Float_t w_2) {

                            if(Gen && Reco){T[i].Fill((double)reco, (double)gen, w_1 * w_2);}

                            else if (Gen && !Reco){T[i].Miss((double)gen);}

                            else if (!Gen && Reco){T[i].Fake((double)reco);}

                        },{Form("gen_%s", vars[i]), Form("reco_%s", vars[i]), "IsGen", "IsReco", "genWeight", "L1PreFiringWeight_Nom"});
    }
    
    //Output file

    TFile output(outname, "RECREATE");

    //Final histogram
    auto h_mc = new_df.Histo1D({"M_invMC", "DiMuon Mass MC", 100, 70, 110}, "gen_mass");
    h_mc->Write("InvMass_MC");
    auto h_reco = new_df.Histo1D({"M_inv_reco", "DiMuon Mass reco", 100, 70, 110}, "reco_mass");
    h_reco->Write("InvMass_reco");

    //Visualization of response matrix and reconstruction efficiency

    for (size_t i = 0; i < vars.size(); i++){

        T[i].Write(Form("Response_%s", vars[i]));

        TH2D* M = (TH2D*) T[i].Hresponse();

        TH1D* h_true = (TH1D*) T[i].Htruth();

        TH1D* h_matched = M->ProjectionY();

        TH1D* h_eff = (TH1D*) h_matched->Clone();

        h_eff->Divide(h_matched, h_true, 1.0, 1.0, "B");

        h_eff->SetLineColor(kBlue);
        h_eff->SetMarkerColor(kBlue);
        h_eff->SetMarkerStyle(20);

        h_eff->Write(Form("Efficiency for %s", vars[i]));

    }

    //Validation for response matrix

    for (size_t i = 0; i < vars.size(); i++){

        if (vars.size() != (bounds.size())){

            std::cerr << "Mismatch : bounds is long " << bounds.size() << " while vars " << vars.size() << std::endl;

        }

        TCanvas* c = new TCanvas(vars[i], vars[i], 800, 600);

        auto h_true_ptr = new_df.Histo1D({Form("%s_MC", vars[i]), Form("%s_MC", vars[i]), 100, bounds[i].first, bounds[i].second}, Form("gen_%s", vars[i]));

        auto h_obt_ptr = new_df.Histo1D({Form("%s_Reco", vars[i]), Form("%s_Reco", vars[i]), 100, bounds[i].first, bounds[i].second}, Form("reco_%s", vars[i]));

        TH1D h_true = *h_true_ptr;
        TH1D h_obt = *h_obt_ptr;

        T[i].UseOverflow();

        RooUnfoldBayes unfold(&T[i], &h_obt, 20);

        TH1D* hUnfold = (TH1D*) unfold.Hunfold();

        hUnfold->SetLineColor(kBlue);
        hUnfold->SetMarkerColor(kBlue);
        hUnfold->SetMarkerStyle(20);

        hUnfold->Draw("P E");

        h_true.Draw("SAME");

        c->Update();

        c->Write(Form("%s_unfolded", vars[i]));

        h_true.Scale(1.0 / h_true.Integral());

        h_true.Write(Form("%s_true", vars[i]));
        
    }

    //Closing output file

    output.Close();

    //Ending Chrono counter e elapsed time;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    //Elapsed time printing

    if (!mute) {std::cout << "Tempo di esecuzione: " << elapsed.count() << " secondi." << std::endl;}

    if (!mute) {std::cout << "Response terminato" << std::endl;}

    //Deactivating batch mode

    gROOT->SetBatch(kFALSE);
}

/**
 * @brief CrossSection wrapper that applies Bayesian unfolding (RooUnfoldBayes) to three reconstructed histograms
 *        using a precomputed response matrix.
 * @param folder_name : folder path, required to contain at least a proper .root file;
 * @param n_iter Number of iterations for the Bayesian unfolding algorithm.
 *               Typical values range from 2 to 10; higher values may introduce fluctuations.
 * @param rpath Path to the ROOT file containing the response matrix (generated with Response.cpp).
 *              The matrix must be compatible in binning and observable definition.
 * @param outname Name of the output ROOT file where unfolded histograms will be stored.
 *
 * @return std::vector<TH1D*> Vector containing the three unfolded histograms, in the following order:
 *         [0] variable Pt, [1] variable Phi_eta, [2] variable y (see CrossSection for definitions).
 *
 * @note Requires ROOT framework and RooUnfold package.
 *
 * @warning Input file must have expected structure, like NanoAD CMS OpenData.
 */

void Unfolded(const char* folder_name,
              int n_iter,
              const char* rpath = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Response.root", 
              const char* outname = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Unfolded.root"){

    //Matrices file

    TFile* Rfile = TFile::Open(rpath);

    if (Rfile->IsZombie()){

        std::cerr << "Response file not opened" << std::endl; 

    }

    //Response matrix vector

    std::vector<RooUnfoldResponse> T(vars.size());

    //Filling matrices

    for(size_t i = 0; i < vars.size(); i++){

        RooUnfoldResponse* M = (RooUnfoldResponse*) Rfile->Get(Form("Response_%s", vars[i]));

        if (!M) {

            std::cerr << "Errore: Response_" << vars[i] << " non trovato!" << std::endl;
            continue;

        }

        T[i] = *M;

    }

    //Closing file

    Rfile->Close();

    //Executing CrossSection

    std::vector<TH1D> h_v = CrossSection(folder_name, "Repo/outFiles/NotUnfolded.root", true, true);

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

        T[i].UseOverflow();

        RooUnfoldBayes unfold(&T[i], h, n_iter);

        TH1D* hUnfold = (TH1D*) unfold.Hunfold();

        if (!hUnfold) {
            std::cerr << "Unfold fallito per i=" << i << std::endl;
            continue;
        }

        h_u[i] = *hUnfold;

        hUnfold->Write(Form("%s_unfolded", vars[i]));

    }

    output->Close();

    std::cout << "Saved correctly" << std::endl;

    std::cout << "Unfolding terminato" << std::endl;

}

/**
 * @brief Compares reconstructed (not unfolded) and unfolded distributions by overlaying the corresponding histograms. 
 *        Also makes a graphs of normalized measured and simulated distributions.
 *        Mainly intended to make stylish graphs and for diagnostic of previous outputs.
 *
 * @param f1 : path to not unfolded histograms, assumed to be output of CrossSection;
 * @param f2 : path to unfolded histograms, assumed to be output of Unfolded;
 * @param f3 : path to MC histograms, assumed to be output of Response
 * @param outname Name of the output ROOT file where comparison plots are saved.
 *
 * @return void. The function produces comparison plots (TCanvas objects)
 *         stored in the output file and accessible via TBrowser.
 *
 * @note Requires ROOT framework and RooUnfold package.
 *
 * @warning Input file must have expected structure, like NanoAD CMS OpenData.
 */

void Comp(const char* f1,
          const char* f2,
          const char* f3,  
          const char* outname = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Comparison.root",
          bool mute = false){

    //Results not-unfolded and unfolded, and MC

    TFile* NU = TFile::Open(f1);

    if (NU->IsZombie()) {std::cerr << "NotUnfolded not opened" << std::endl;}

    std::vector<TH1D> h_nu(vars.size());

    TFile* U = TFile::Open(f2);

    if (U->IsZombie()) {std::cerr << "Unfolded not opened" << std::endl;}

    std::vector<TH1D> h_u(vars.size());

    TFile* MC = TFile::Open(f3);

    if (MC->IsZombie()) {std::cerr << "Montecarlo not opened" << std::endl;}

    std::vector<TH1D> h_mc(vars.size());

    for (size_t i = 0; i < vars.size(); i++){

        TH1D* h = (TH1D*) NU->Get(vars[i]); 

        h_nu[i] = *h;

        TH1D* hu = (TH1D*) U->Get(Form("%s_unfolded", vars[i]));

        h_u[i] = *hu;

        TH1D* hMC = (TH1D*) MC->Get(Form("%s_true", vars[i]));

        h_mc[i] = *hMC;

    }

    //Outfile

    TFile* output = new TFile(outname, "RECREATE");

    if (h_nu.size() == h_u.size()){

        for (size_t i = 0; i < h_u.size(); i++){

            TCanvas* c = new TCanvas(vars[i], vars[i], 800, 600);

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

            TH1D h_D = h_u[i]; 

            h_u[i].Scale(1./(L * h_u[i].GetXaxis()->GetBinWidth(1)));

            h_u[i].Write(Form("CrossSection_%s", vars[i]));

            TCanvas* c_MC = new TCanvas(Form("%s_MC", vars[i]), Form("%s_MC", vars[i]), 800, 600);

            TLegend* legend_MC = new TLegend(0.7, 0.7, 0.9, 0.9);

            h_mc[i].SetStats(0);

            h_mc[i].Draw();

            legend_MC->AddEntry(&h_mc[i], "Montecarlo", "lep");

            h_D.Scale(1./(s*L));

            legend_MC->AddEntry(&h_D, "Data", "lep");

            h_D.Draw("SAME");

            legend_MC->Draw("SAME");

            c_MC->Update();

            c_MC->Write();
        }

    }

    output->Close();

}

