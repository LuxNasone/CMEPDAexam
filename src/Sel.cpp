//Standard include

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

//ROOT include

#include <TSystem.h> 
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>

//Useful functions:

//Invariant mass calculator

ROOT::RVec<Float_t> Minv_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p){
    ROOT::RVec<Float_t> Minv;
    for (size_t i = 0; i < p.size(); i++){for (size_t j = i + 1; j < p.size(); ++j){Minv.push_back((p[i] + p[j]).M());}}
    return Minv;
}

//Selection at generator level (used in MCSel)

ROOT::RVec<ROOT::Math::PtEtaPhiMVector> GenSel(const UInt_t &n, const ROOT::RVec<int> &id, const ROOT::RVec<int> &mid, const ROOT::RVecF &pt, const ROOT::RVecF &eta, const ROOT::RVecF &phi, const ROOT::RVecF &mass){

    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> vecs;
    vecs.reserve(n);

    for (UInt_t i = 0; i < n; i++){
                                                        
        if ((id[i] == 13 || id[i] == -13) && mid[i] >= 0 && id[mid[i]] == 23){

            ROOT::Math::PtEtaPhiMVector p(pt[i],eta[i],phi[i], mass[i]);

            vecs.push_back(p);

            }
        }
        return vecs;
}

//Selection applied in the article

ROOT::RVec<ROOT::Math::PtEtaPhiMVector> Reco(const UInt_t &n, const ROOT::RVec<int> &charge, const ROOT::RVec<float> &Iso,const ROOT::RVec<float> &pt, const ROOT::RVec<float> &eta, const ROOT::RVec<float> phi,const ROOT::RVec<float> &mass){

    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> vecs;

    if (n == 2){
        if (charge[0] + charge[1] == 0 && Iso[0] < 0.15 && Iso[1] < 0.15 && pt[0] > 25 && pt[1] > 25 && abs(eta[0]) < 2.4 && abs(eta[1]) < 2.4){

            vecs.push_back(ROOT::Math::PtEtaPhiMVector(pt[0],eta[0],phi[0], mass[0]));

            vecs.push_back(ROOT::Math::PtEtaPhiMVector(pt[1],eta[1], phi[1], mass[1]));

        }
    }

    return vecs;

}

//Selection on invariant mass

bool Minv_Range(const ROOT::RVec<Float_t>& masses) {
    for (auto m : masses) {if (std::abs(m - 91.1817) < 15) {return true;}}
    return false;
}

//Macro to reconstruct Z^0 peak and obtain histograms of pt, eta and phi;

void DataSel(const char* fname, std::string outname = "Repo/outFiles/OutDS.root", bool MT = true){

    //Necessary imports for 4-vectors and other utilities, compile macro with + at the end;

    gSystem->Load("libPhysics"); 
    gSystem->Load("libMathCore");
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

    */

    auto new_df = df.Filter("nMuon == 2 && "
                              "Muon_charge[0] + Muon_charge[1] == 0 &&"
                              "Muon_pfRelIso03_all[0] < 0.15 && Muon_pfRelIso03_all[1] < 0.15 &&"
                              "Muon_pt[0] > 25 && Muon_pt[1] > 25 &&"
                              "abs(Muon_eta[0]) < 2.4 && abs(Muon_eta[1]) < 2.4")
                      .Define("Dimuon_mass", [](const ROOT::RVecF &pt, const ROOT::RVecF &eta, const ROOT::RVecF &phi, const ROOT::RVecF &mass){
                               ROOT::Math::PtEtaPhiMVector p_1(pt[0],eta[0],phi[0], mass[0]), p_2(pt[1],eta[1], phi[1], mass[1]);
                               return (Float_t)(p_1 + p_2).M();}, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"})
                      .Filter("abs(Dimuon_mass - 91.1817) < 15");

    std::vector<std::string> vars = {"pt", "eta", "phi"};

    for (int i = 0; i < 2; i++){

        for (auto v : vars){

            new_df = new_df.Define(Form("Muon%d_%s", i, v.c_str()), Form("Muon_%s[%d]", v.c_str(), i));

            }
        }

    //Output file

    TFile output(outname.c_str(), "RECREATE");

    //Histogram for Invariant Mass;

    auto h_dmm = new_df.Histo1D({"M_inv", "DiMuon_mass", 100, 70, 110}, "Dimuon_mass");
    h_dmm->Draw();
    h_dmm->Write("InvMass");


    //Bounds for drawing histograms;
    std::vector<std::pair<Float_t, Float_t>> bounds = {{15, 200}, {-3, 3}, {-4, 4}};

    //Loop for diMuon kinematic variables histograms;

    for (int i = 0; i < 3; i++){

        auto h_1 = new_df.Histo1D(ROOT::RDF::TH1DModel(vars[i].c_str(), (vars[i] + to_string(1)).c_str(), 100, bounds[i].first, bounds[i].second), 
                                   Form("Muon%d_%s", 0, vars[i].c_str()));
        auto h_2 = new_df.Histo1D(ROOT::RDF::TH1DModel(vars[i].c_str(), (vars[i] + to_string(2)).c_str(), 100, bounds[i].first, bounds[i].second), 
                                   Form("Muon%d_%s", 1, vars[i].c_str()));


        h_1->Draw();
        h_2->Draw("SAME");

        h_1->Write((vars[i] + to_string(1)).c_str());
        h_2->Write((vars[i] + to_string(2)).c_str());
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
    gROOT->SetBatch(kTRUE);

    //Chrono counter;

    auto start = std::chrono::high_resolution_clock::now();

    //Optional: activate multithreading

    if (MT){ROOT::EnableImplicitMT();}

    //Initializing DataFrame, fname must be to a .root file;

    ROOT::RDataFrame df("Events", fname);
    
    /* Reconstructing Z^0 peak by MC tagging at generator level:

    Loop on GenPart_pdgId, need at least a +13 and -13 con GenPart_pdgId(GenPart_genPartIdxMother) == 23 

    We then compute invariant mass for the pair

    */

    auto new_df = df.Define("Muon_From_Z", GenSel, {"nGenPart", "GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass"})
                    .Define("WR", Reco, {"nMuon", "Muon_charge", "Muon_pfRelIso03_all","Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"})
                    .Define("DiMuonMass", Minv_calculator, {"Muon_From_Z"})
                    .Filter(Minv_Range, {"DiMuonMass"})
                    .Define("WRMass", Minv_calculator, {"WR"})
                    .Filter(Minv_Range, {"WRMass"});

    std::vector<std::string> vars = {"pt", "eta", "phi"};

    for (int i = 0; i < 2; i++){

        for (auto v : vars){

            new_df = new_df.Define(Form("Muon%d_%s", i, v.c_str()), Form("Muon_From_Z[%d].%s()", i, v.c_str()));
            new_df = new_df.Define(Form("Muon%d%sWR", i, v.c_str()), Form("WR[%d].%s()", i, v.c_str()));

            }
        }
                    

    //Output file

    TFile output(outname.c_str(), "RECREATE");

    //Final histogram
    auto h_mc = new_df.Histo1D({"M_invMC", "DiMuon Mass MC", 100, 70, 110}, "DiMuonMass");
    h_mc->Draw();
    h_mc->Write("InvMass_MC");

    //Bounds for drawing histograms;
    std::vector<std::pair<Float_t, Float_t>> bounds = {{15, 200}, {-3, 3}, {-4, 4}};

    //Loop for diMuon kinematic variables histograms;

    for (int i = 0; i < 3; i++){

        auto h_1 = new_df.Histo1D(ROOT::RDF::TH1DModel(vars[i].c_str(), (vars[i] + to_string(1)).c_str(), 100, bounds[i].first, bounds[i].second), 
                                  Form("Muon%d_%s", 0, vars[i].c_str()));
        auto h_2 = new_df.Histo1D(ROOT::RDF::TH1DModel(vars[i].c_str(), (vars[i] + to_string(2)).c_str(), 100, bounds[i].first, bounds[i].second), 
                                 Form("Muon%d_%s", 1, vars[i].c_str()));

        h_1->Draw();
        h_2->Draw("SAME");

        h_1->Write((vars[i] + to_string(1) + "_MC").c_str());
        h_2->Write((vars[i] + to_string(2) + "_MC").c_str());

    }

    //Response matrix

    auto T = new_df.Histo2D({"hT", "Response Matrix", 100, 70, 110, 100, 70, 110}, "WRMass", "DiMuonMass"); 
    T->Draw();
    T->Write("InvMass_Response");


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