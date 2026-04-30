#include <cmath>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>

int n_b = 30; 

std::vector<const char*> vars = {"pt", "phi_eta", "y"};

std::vector<const char*> xlabels = {"p^{Z}_{T}[GeV]", "#phi^{*}_{#eta}", "|y^{Z}|"};

std::vector<const char*> ylabels = {"d#sigma / dp^{Z}_{T}[pb/GeV]", "d#sigma / d#phi^{*}_{#eta} [pb]", "d#sigma / dy^{Z} [pb]"};

std::vector<const char*> titles = {"Transverse momentum", "Optimized angle", "Rapidity abs."};
 
std::vector<std::pair<Float_t, Float_t>> bounds = {{1, 100}, {0.01, 3}, {0.01, 2.5}};

std::vector<std::pair<Float_t, Float_t>> range = {{0, 8.5e5}, {0, 1.2e6}, {0, 1.2e5}};

double L = 8740.119304;

double s = 1952;

Double_t Minv_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p){

    if (p.size()>= 2){return (p[0] + p[1]).M();}

    return 0;
}

Double_t Pt_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p){

    if (p.size()>= 2){return (p[0] + p[1]).Pt();}

    return 0;
}

Double_t y_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p){

    if (p.size()>= 2){return abs((p[0] + p[1]).Rapidity());}

    return 0;

}

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

bool IsReco(const UInt_t &n, const ROOT::RVec<int> &charge, const ROOT::RVec<float> &Iso,const ROOT::RVec<float> &pt, const ROOT::RVec<float> &eta){
    
    bool IsReco = false;

    if (n == 2){
        if (charge[0] + charge[1] == 0 && Iso[0] < 0.15 && Iso[1] < 0.15 && pt[0] > 25 && pt[1] > 25 && abs(eta[0]) < 2.4 && abs(eta[1]) < 2.4){

            IsReco = true;

        }
    }

    return IsReco;
}

ROOT::RVec<ROOT::Math::PtEtaPhiMVector> GenSel(const UInt_t &n, const ROOT::RVec<int> &id, const ROOT::RVec<int> &mid, const ROOT::RVec<float> &pt, const ROOT::RVec<float> &eta, const ROOT::RVec<float> &phi, const ROOT::RVec<float> &mass){

    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> vecs;

    for (UInt_t i = 0; i < n; i++){
        if (!(id[i] == 13 || id[i] == -13)) continue;
        if (mid[i] < 0 || id[mid[i]] != 23) continue;

        for (UInt_t j = i+1; j < n; j++){
            if (!(id[j] == 13 || id[j] == -13)) continue;
            if (mid[j] < 0 || id[mid[j]] != 23) continue;

            if (id[i] == -id[j]) {
                vecs.emplace_back(pt[i], eta[i], phi[i], mass[i]);
                vecs.emplace_back(pt[j], eta[j], phi[j], mass[j]);
                return vecs; 
            }
        }
    }

    return vecs;
}

ROOT::RVec<ROOT::Math::PtEtaPhiMVector> Reco(const UInt_t &n, const ROOT::RVec<int> &charge, const ROOT::RVec<float> &Iso,const ROOT::RVec<float> &pt, const ROOT::RVec<float> &eta, const ROOT::RVec<float> &phi,const ROOT::RVec<float> &mass){

    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> vecs;

    if (n == 2){
        if (charge[0] + charge[1] == 0 && Iso[0] < 0.15 && Iso[1] < 0.15 && pt[0] > 25 && pt[1] > 25 && abs(eta[0]) < 2.4 && abs(eta[1]) < 2.4){

            vecs.push_back(ROOT::Math::PtEtaPhiMVector(pt[0],eta[0],phi[0], mass[0]));

            vecs.push_back(ROOT::Math::PtEtaPhiMVector(pt[1],eta[1], phi[1], mass[1]));

        }
    }

    return vecs;

}

bool Minv_Range(Double_t m) {
    if (std::abs(m - 91.1817) < 15) {return true;}
    return false;
}
