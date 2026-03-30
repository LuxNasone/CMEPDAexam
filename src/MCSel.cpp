#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

#include <TSystem.h> 
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>

//Reconstruction at generator level of Z^0 decays

void MCSel(const char* fname, std::string outname, bool MT = true){

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

    auto new_df = df.Define("Muon_From_Z", [] (const UInt_t &n, const ROOT::RVec<int> &id, const ROOT::RVec<int> &mid, const ROOT::RVecF &pt,
                                                const ROOT::RVecF &eta, const ROOT::RVecF &phi, const ROOT::RVecF &mass){

                                                    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> vecs;
                                                    vecs.reserve(n);

                                                    for (UInt_t i = 0; i < n; i++){
                                                        
                                                        if ((id[i] == 13 || id[i] == -13) && mid[i] >= 0 && id[mid[i]] == 23){

                                                            ROOT::Math::PtEtaPhiMVector p(pt[i],eta[i],phi[i], mass[i]);

                                                            vecs.push_back(p);

                                                        }
                                                    }
                                                    return vecs;
                                                }, {"nGenPart", "GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass"})
                        .Define("DiMuonMass", [](const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p){
                                                    ROOT::RVec<Float_t> Minv;
                                                    for (size_t i = 0; i < p.size(); i++){
                                                        for (size_t j = i + 1; j < p.size(); ++j){Minv.push_back((p[i] + p[j]).M());}}
                                                    return Minv;
                        }, {"Muon_From_Z"})
                        .Filter([](const ROOT::RVec<Float_t>& masses) {
                                    for (auto m : masses) {
                                        if (std::abs(m - 91.1817) < 15) {return true;}
                                        }
                                    return false;
                                    }, {"DiMuonMass"})
                        .Define("Muon0_pt",  "Muon_From_Z[0].Pt()")
                        .Define("Muon1_pt",  "Muon_From_Z[1].Pt()")
                        .Define("Muon0_eta", "Muon_From_Z[0].Eta()")
                        .Define("Muon1_eta", "Muon_From_Z[1].Eta()")
                        .Define("Muon0_phi", "Muon_From_Z[0].Phi()")
                        .Define("Muon1_phi", "Muon_From_Z[1].Phi()");

    //Final histogram
    auto h_mc = new_df.Histo1D({"M_invMC", "DiMuon Mass MC", 100, 70, 110}, "DiMuonMass");

    TCanvas c("M_invMC_canvas", "DiMuon Mass MC", 800, 600);
    h_mc->Draw();
    c.SaveAs((outname + "MCMass.png").c_str());

    //Variables necessary for drawing histograms;

    std::vector<std::string> vars = {"pt", "eta", "phi"};
    std::vector<std::pair<Float_t, Float_t>> bounds = {{15, 200},
                                              {-3, 3},
                                              {-4, 4}};

    //Loop for diMuon kinematic variables histograms;

    for (int i = 0; i < 3; i++){

        TCanvas *c = new TCanvas((vars[i] + "canvas").c_str(), vars[i].c_str(), 800, 600);

        auto h_1 = new_df.Histo1D(ROOT::RDF::TH1DModel(vars[i].c_str(), (vars[i] + to_string(1)).c_str(), 100 ,bounds[i].first, bounds[i].second), 
                                  Form("Muon%d_%s", 0, vars[i].c_str()));
        auto h_2 = new_df.Histo1D(ROOT::RDF::TH1DModel(vars[i].c_str(), (vars[i] + to_string(2)).c_str(), 100 , bounds[i].first, bounds[i].second), 
                                 Form("Muon%d_%s", 1, vars[i].c_str()));

        h_1->SetLineColor(kBlue);
        h_2->SetLineColor(kRed);

        h_1->Draw();
        h_2->Draw("SAME");

        c->Update();
        c->SaveAs((outname + vars[i] + "MC.png").c_str());

        delete c;
    }

    //Ending Chrono counter e elapsed time;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    //Elapsed time printing

    std::cout << "Tempo di esecuzione: " << elapsed.count() << " secondi." << std::endl;


}