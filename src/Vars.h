#ifndef VARS_H
#define VARS_H

#include <Rtypes.h>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>

Double_t Minv_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p);

Double_t Pt_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p);

Double_t y_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p);

Double_t phi_eta_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p);

ROOT::RVec<ROOT::Math::PtEtaPhiMVector> GenSel(const UInt_t &n, const ROOT::RVec<int> &id, const ROOT::RVec<int> &mid, const ROOT::RVec<float> &pt, const ROOT::RVec<float> &eta, const ROOT::RVec<float> &phi, const ROOT::RVec<float> &mass);

ROOT::RVec<ROOT::Math::PtEtaPhiMVector> Reco(const UInt_t &n, const ROOT::RVec<int> &charge, const ROOT::RVec<float> &Iso,const ROOT::RVec<float> &pt, const ROOT::RVec<float> &eta, const ROOT::RVec<float> &phi,const ROOT::RVec<float> &mass);

bool Minv_Range(Double_t m);

bool IsReco(const UInt_t &n, const ROOT::RVec<int> &charge, const ROOT::RVec<float> &Iso,const ROOT::RVec<float> &pt, const ROOT::RVec<float> &eta);


bool IsTrue(const UInt_t &n_g,const ROOT::RVec<int> &id, const ROOT::RVec<int> &mid);

#endif