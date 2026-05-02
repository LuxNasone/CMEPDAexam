#include <iostream>
#include <fstream>
#include <cmath>
#include <numbers>
#include <chrono>

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

#include "../include/Utils.h"

std::vector<TH1D> NotUnfolded(const char* folder_name,
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

    ROOT::RDataFrame df("Events", Form("%s/*.root",folder_name));

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
                              "abs(Muon_eta[0]) < 2.4 && abs(Muon_eta[1]) < 2.4 &&"
                              "abs(Muon_mass[0] - 0.1057) < 2.5e-5  && abs(Muon_mass[1] - 0.1057) < 2.5e-5")
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

        auto h_tmp = new_df.Histo1D(ROOT::RDF::TH1DModel(vars[i], titles[i], n_b, bounds[i].first, bounds[i].second), Form("%s", vars[i]));

        TH1D h = h_tmp.GetValue(); 

        h_v[i] = h;

    }

    //Saving on file

    TFile output(outname, "RECREATE");

    //Histogram for Invariant Mass;

    auto h_m = new_df.Histo1D({"M_inv", "Z^{0} mass", n_b, 70, 110}, "mass");

    //Styling

    TCanvas* c_m = new TCanvas("M_inv", "Z^{0} mass", 800, 600);

    h_m->GetXaxis()->SetTitle("M_{inv} [GeV]");

    h_m->GetYaxis()->SetTitle("Counts [pure]");

    h_m->SetLineColor(kBlue);

    h_m->SetMarkerColor(kBlue);

    h_m->SetMarkerStyle(20);

    h_m->SetStats(0);
    
    h_m->Draw("P E");

    c_m->Update();

    h_m->Write("InvMass");

    c_m->Write("InvMassCanvas");

    //Writing for other histograms and styling

    for (size_t i = 0; i < vars.size(); i++){

        TCanvas* c_v = new TCanvas(vars[i], vars[i], 800, 600);

        h_v[i].GetXaxis()->SetTitle(xlabels[i]);

        h_v[i].GetYaxis()->SetTitle("Counts [pure]");

        h_v[i].SetLineColor(kBlue);

        h_v[i].SetMarkerColor(kBlue);

        h_v[i].SetMarkerStyle(20);

        h_v[i].SetStats(0);
        
        h_v[i].Draw("P E");

        h_v[i].Write(vars[i]);

        c_v->Write(Form("%s_Canvas", vars[i]));

        delete c_v;

    }

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

    ROOT::RDataFrame df("Events", Form("%s/*.root",folder_name));

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
                    .Define("IsReco", IsReco, {"nMuon", "Muon_charge", "Muon_pfRelIso03_all","Muon_pt", "Muon_eta", "Muon_mass"})
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

        new_df.Foreach([&](Double_t gen, Double_t reco, bool Gen, bool Reco, Float_t w_1) {

                            if(Gen && Reco){T[i].Fill((double)reco, (double)gen, w_1);}

                            else if (Gen && !Reco){T[i].Miss((double)gen);}

                            else if (!Gen && Reco){T[i].Fake((double)reco);}

                        },{Form("gen_%s", vars[i]), Form("reco_%s", vars[i]), "IsGen", "IsReco", "genWeight"});
    }
    
    //Output file

    TFile output(outname, "RECREATE");

    //Final histogram
    auto h_mc = new_df.Histo1D({"M_invMC", "Z^{0} mass (MC)", n_b, 70, 110}, "gen_mass");

    TCanvas* c_mc = new TCanvas("M_invMC", "Z^{0} mass (MC)", 800, 600);

    h_mc->GetXaxis()->SetTitle("M_{inv} [GeV]");

    h_mc->GetYaxis()->SetTitle("Counts [pure]");

    h_mc->SetLineColor(kBlue);

    h_mc->SetMarkerColor(kBlue);

    h_mc->SetMarkerStyle(20);

    h_mc->SetMarkerSize(1.2);

    h_mc->SetStats(0);
    
    h_mc->Draw("P E");

    c_mc->Update();

    h_mc->Write("InvMass_MC");

    c_mc->Write("InvMass_MCcanvas");

    auto h_reco = new_df.Histo1D({"M_inv_reco", "Z^{0} mass (Reco)", n_b, 70, 110}, "reco_mass");

    TCanvas* c_reco = new TCanvas("M_invReco", "Z^{0} mass (Reco)", 800, 600);

    h_reco->GetXaxis()->SetTitle("M_{inv} [GeV]");

    h_reco->GetYaxis()->SetTitle("Counts [pure]");

    h_reco->SetLineColor(kBlue);

    h_reco->SetMarkerColor(kBlue);

    h_reco->SetMarkerStyle(20);

    h_reco->SetMarkerSize(1.2);

    h_reco->SetStats(0);
    
    h_reco->Draw("P E");

    c_reco->Update();

    h_reco->Write("InvMass_reco");

    c_reco->Write("InvMass_MCcanvas");

    //Visualization of response matrix and reconstruction efficiency

    for (size_t i = 0; i < vars.size(); i++){

        T[i].Write(Form("Response_%s", vars[i]));

        TCanvas* c_T = new TCanvas(Form("Response : %s", titles[i]), Form("Response : %s", titles[i]), 800, 600);

        TH2D* M = (TH2D*) T[i].Hresponse();

        M->SetTitle(Form("Response : %s", titles[i]));

        M->GetXaxis()->SetTitle(Form("%s (reconstructed)", xlabels[i]));

        M->GetYaxis()->SetTitle(Form("%s (generated)", xlabels[i]));

        M->SetStats(0);

        M->Draw("COLZ");

        c_T->Update();

        M->Write(Form("Response matrix for %s", titles[i]));

        c_T->Write(Form("Response matrix for %s canvas", titles[i]));

        delete c_T;

        TCanvas* c_eff = new TCanvas(Form("Efficiency for %s", titles[i]), Form("Efficiency for %s", titles[i]), 800, 600);

        TH1D* h_true = (TH1D*) T[i].Htruth();

        TH1D* h_matched = M->ProjectionY();

        TH1D* h_eff = (TH1D*) h_matched->Clone();

        h_eff->Divide(h_matched, h_true, 1.0, 1.0, "B");

        h_eff->SetTitle(Form("Reconstruction efficiency for %s", titles[i]));

        h_eff->GetXaxis()->SetTitle(xlabels[i]);

        h_eff->GetYaxis()->SetTitle("Efficency [pure]");

        h_eff->SetLineColor(kBlue);

        h_eff->SetMarkerColor(kBlue);

        h_eff->SetMarkerStyle(20);

        h_eff->SetStats(0);

        h_eff->SetMarkerSize(1.2);

        h_eff->Draw();

        c_eff->Update();

        h_eff->Write(Form("Efficiency for %s", vars[i]));

        c_eff->Write(Form("Efficiency for %s (canvas)", vars[i]));

        delete c_eff;

    }

    //Validation for response matrix

    for (size_t i = 0; i < vars.size(); i++){

        if (vars.size() != (bounds.size())){

            std::cerr << "Mismatch : bounds is long " << bounds.size() << " while vars " << vars.size() << std::endl;

        }

        TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);

        auto h_true_ptr = new_df.Histo1D({Form("%s (MC)", titles[i]), Form("%s (MC)", titles[i]), n_b, bounds[i].first, bounds[i].second}, Form("gen_%s", vars[i]));

        auto h_obt_ptr = new_df.Histo1D({Form("%s (Reco)", titles[i]), Form("%s (Reco)", titles[i]), n_b, bounds[i].first, bounds[i].second}, Form("reco_%s", vars[i]));

        TH1D h_true = *h_true_ptr;

        TCanvas* c_true = new TCanvas(Form("%s (MC)", titles[i]), Form("%s (MC)", titles[i]), 800, 600);

        h_true.GetXaxis()->SetTitle(xlabels[i]);

        h_true.GetYaxis()->SetTitle("Counts [pure]");

        h_true.SetLineColor(kBlue);

        h_true.SetMarkerColor(kBlue);

        h_true.SetMarkerStyle(20);

        h_true.SetStats(0);

        h_true.Draw("P E");

        c_true->Update();

        c_true->Write(Form("%s true (canvas)", vars[i]));

        leg->AddEntry(&h_true, "True", "lep");

        TCanvas* c_obt = new TCanvas(Form("%s (Reco)", titles[i]), Form("%s (Reco)", titles[i]), 800, 600);
        
        TH1D h_obt = *h_obt_ptr;

        h_obt.GetXaxis()->SetTitle(xlabels[i]);

        h_obt.GetYaxis()->SetTitle("Counts [pure]");

        h_obt.SetLineColor(kBlue);

        h_obt.SetMarkerColor(kBlue);

        h_obt.SetMarkerStyle(20);

        h_obt.SetStats(0);

        h_obt.Draw("P E");

        c_obt->Update();

        c_obt->Write(Form("%s obt (canvas)", vars[i]));

        T[i].UseOverflow();

        TCanvas* c = new TCanvas(vars[i], vars[i], 800, 600);

        RooUnfoldBayes unfold(&T[i], &h_obt, 5);

        TH1D* hUnfold = (TH1D*) unfold.Hunfold();

        hUnfold->GetXaxis()->SetTitle(xlabels[i]);

        hUnfold->GetYaxis()->SetTitle("Counts [pure]");

        hUnfold->SetTitle(Form("Comparison generated unfolded (%s)", titles[i]));

        hUnfold->SetLineColor(kRed);

        hUnfold->SetMarkerColor(kRed);

        hUnfold->SetMarkerStyle(20);

        hUnfold->SetStats(0);

        leg->AddEntry(hUnfold, "Unfolded", "lep");

        hUnfold->Draw("P E");

        h_true.Draw("P E SAME");

        leg->Draw("SAME");

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

    //Executing NotUnfolded

    std::vector<TH1D> h_v = NotUnfolded(folder_name, "Repo/outFiles/NotUnfolded.root", true, true);

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

void Comparison(const char* f1,
          const char* f2,
          const char* f3,  
          const char* outname = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Comparison.root",
          bool mute = false){

    gROOT->SetBatch(kTRUE);

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

            h_nu[i].Draw("P E");

            h_u[i].SetLineColor(kRed);

            h_u[i].SetMarkerColor(kRed);

            h_u[i].SetMarkerStyle(20);

            h_u[i].SetStats(0);

            legend->AddEntry(&h_u[i], "Unfolded", "lep");

            h_u[i].Draw("P E SAME");

            legend->Draw();

            c->Update();

            c->Write();

            TH1D h_D = h_nu[i]; 

            TCanvas* c_CS = new TCanvas(Form("CrossSection_%s", vars[i]), Form("CrossSection_%s", vars[i]), 800, 600);

            h_u[i].Scale(1./(L * h_u[i].GetXaxis()->GetBinWidth(1)));

            h_u[i].GetXaxis()->SetTitle(xlabels[i]);

            h_u[i].GetYaxis()->SetTitle(ylabels[i]);

            h_u[i].SetTitle(Form("Differential Cross Section (%s)", titles[i]));

            h_u[i].Draw("P E");

            c_CS->Update();

            c_CS->Write(Form("CrossSection_%s", vars[i]));

            TCanvas* c_MC = new TCanvas(Form("%s_MC", vars[i]), Form("%s_MC", vars[i]), 800, 600);

            TLegend* legend_MC = new TLegend(0.7, 0.7, 0.9, 0.9);

            h_mc[i].GetXaxis()->SetTitle(xlabels[i]);

            h_mc[i].GetYaxis()->SetTitle("Norm. counts [pure]");

            h_mc[i].SetTitle(Form("Comparison between MC and measured (%s)", titles[i]));

            h_mc[i].SetLineColor(kRed);

            h_mc[i].SetMarkerColor(kRed);

            h_mc[i].SetMarkerStyle(20);

            h_mc[i].SetStats(0);

            h_mc[i].Draw();

            legend_MC->AddEntry(&h_mc[i], "Montecarlo", "lep");

            h_D.Scale(1./h_D.Integral());

            legend_MC->AddEntry(&h_D, "Data", "lep");

            h_D.Draw("SAME");

            legend_MC->Draw("SAME");

            c_MC->Update();

            c_MC->Write();
        }

    }

    output->Close();

    gROOT->SetBatch(kFALSE);

}
