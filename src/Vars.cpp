#include <cmath>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>

//Useful functions:

//Invariant mass calculator

Double_t Minv_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p){

    if (p.size()>= 2){return (p[0] + p[1]).M();}

    return 0;
}

//Z^0 transverse momentum

Double_t Pt_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p){

    if (p.size()>= 2){return (p[0] + p[1]).Pt();}

    return 0;
}

//Z^0 rapidity

Double_t y_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p){

    if (p.size()>= 2){return abs((p[0] + p[1]).Rapidity());}

    return 0;

}

//Z^0 phi*eta
Double_t phi_eta_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p){

    if (p.size()>= 2){

            double dphi = std::fabs(p[0].Phi() - p[1].Phi());
            if (dphi > M_PI) dphi = 2.0 * M_PI - dphi;

            double deta = p[0].Eta() - p[1].Eta();

            double sinThetaStar = 1.0 / std::cosh(deta / 2.0);
            double tanTerm = std::tan((M_PI - dphi) / 2.0);

            Float_t phi_eta = tanTerm * sinThetaStar;

            return phi_eta;
    }

    return 0;

}

//Function to understand if the event is a true event

bool IsTrue(const UInt_t &n_g,const ROOT::RVec<int> &id, const ROOT::RVec<int> &mid){

    bool pair = false;

    std::vector<int> flag;

    for (UInt_t i = 0; i < n_g; i++){

        if ((id[i] == 13 || id[i] == -13) && mid[i] >= 0 && id[mid[i]] == 23){flag.push_back(id[i]);}

    }

    if (flag.size() >= 2){

        for (size_t i = 0; i < flag.size(); i++){

            for (size_t j = i+1; j < flag.size(); j++){

                if (flag[i] == -flag[j]){pair = true;} 

            }

        }

    }
    
    return pair;

}

//Function to understand if the event is a reconstructed event

bool IsReco(const UInt_t &n, const ROOT::RVec<int> &charge, const ROOT::RVec<float> &Iso,const ROOT::RVec<float> &pt, const ROOT::RVec<float> &eta){
    
    bool IsReco = false;

    if (n == 2){
        if (charge[0] + charge[1] == 0 && Iso[0] < 0.15 && Iso[1] < 0.15 && pt[0] > 25 && pt[1] > 25 && abs(eta[0]) < 2.4 && abs(eta[1]) < 2.4){

            IsReco = true;

        }
    }

    return IsReco;
}

//Selection at generator level (used in MCSel)

ROOT::RVec<ROOT::Math::PtEtaPhiMVector> GenSel(const UInt_t &n, const ROOT::RVec<int> &id, const ROOT::RVec<int> &mid, const ROOT::RVecF &pt, const ROOT::RVecF &eta, const ROOT::RVecF &phi, const ROOT::RVecF &mass){

    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> vecs;
    std::vector<int> flag;
    vecs.reserve(n);

    for (UInt_t i = 0; i < n; i++){
                                                        
        if ((id[i] == 13 || id[i] == -13) && mid[i] >= 0 && id[mid[i]] == 23){

            ROOT::Math::PtEtaPhiMVector p(pt[i],eta[i],phi[i], mass[i]);

            vecs.push_back(p);

            flag.push_back(id[i]);

            }
        }

        bool pair = false;

        if (flag.size() >= 2){

            for (size_t i = 0; i < flag.size(); i++){

                for (size_t j = i+1; j < flag.size(); j++){

                    if (flag[i] == -flag[j]){pair = true;} 

                }
                
                if (pair){break;}

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

bool Minv_Range(Double_t &m) {
    if (std::abs(m - 91.1817) < 15) {return true;}
    return false;
}
