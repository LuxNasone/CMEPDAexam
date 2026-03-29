#include <iostream>
#include <fstream>
#include <cmath>

#include <TSystem.h> 
#include <TFile.h>
#include <TTree.h>
#include <Rtypes.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <ROOT/RDataFrame.hxx>

//Reconstruction at generator level of Z^0 decays

std::vector<std::vector<Float_t>> MCSel(const char* fname, std::string outname, bool MT = true){

    //Necessary for 4 vectors and other utilities, compile with + at the end

    gSystem->Load("libPhysics");
    gSystem->Load("libMathCore");
    gROOT->SetBatch(kTRUE);

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

                                                    for (UInt_t i = 0; i < n; i++){
                                                        
                                                        if ((id[i] == 13 || id[i] == -13) && mid[i] >= 0 && id[mid[i]] == 23){

                                                            ROOT::Math::PtEtaPhiMVector p(pt[i],eta[i],phi[i], mass[i]);

                                                            vecs.push_back(p);

                                                        }
                                                    }
                                                    return vecs;
                                                }, {"nGenPart", "GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass"})
                        .Define("DiMuonMass", [](const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p,  const ROOT::RVec<int> &id, const ROOT::RVec<int> &mid){
                                                    ROOT::RVec<Float_t> Minv;
                                                    for (size_t i = 0; i < p.size(); i++){
                                                        for (size_t j = i + 1; j < p.size(); ++j){
                                                            if (id[i] + id[j] == 0 && mid[i] == mid[j]){
                                                                Minv.push_back((p[i] + p[j]).M());
                                                            }
                                                        }
                                                    }
                                                    return Minv;
                        }, {"Muon_From_Z", "GenPart_pdgId", "GenPart_genPartIdxMother"})
                        .Filter([](const ROOT::RVec<Float_t>& masses) {
                                    for (auto m : masses) {
                                        if (std::abs(m - 91.1817) < 15) {return true;}
                                        }
                                    return false;
                                    }, {"DiMuonMass"});

    //Final histogram
    auto h_mc = new_df.Histo1D({"M_invMC", "DiMuon Mass MC", 100, 70, 110}, "DiMuonMass");

    TCanvas c("M_invMC_canvas", "DiMuon Mass MC", 800, 600);
    h_mc->Draw();
    c.SaveAs((outname + "MCDiMuonMass.png").c_str());

}