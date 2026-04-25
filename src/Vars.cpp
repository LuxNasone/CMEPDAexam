#include <cmath>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>

//Useful functions:

/**
 * @brief Computes the invariant mass of the sum of the first two four-momenta
 *        in the input vector.
 *
 * @param p Vector of four-momenta (e.g. TLorentzVector) representing
 *          reconstructed or generated muons.
 *
 * @return Invariant mass of p[0] + p[1]. Returns 0.0 if the vector size is less than 2.
 *
 */
Double_t Minv_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p){

    if (p.size()>= 2){return (p[0] + p[1]).M();}

    return 0;
}

/**
 * @brief Computes the transverse mmomentum of the sum of the first two four-momenta
 *        in the input vector.
 *
 * @param p Vector of four-momenta (e.g. TLorentzVector) representing
 *          reconstructed or generated muons.
 *
 * @return Transverse momentum of p[0] + p[1]. Returns 0.0 if the vector size is less than 2.
 *
 */

Double_t Pt_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p){

    if (p.size()>= 2){return (p[0] + p[1]).Pt();}

    return 0;
}
/**
 * @brief Computes the rapidity of the sum of the first two four-momenta
 *        in the input vector.
 *
 * @param p Vector of four-momenta (e.g. TLorentzVector) representing
 *          reconstructed or generated muons.
 *
 * @return Rapidity of p[0] + p[1]. Returns 0.0 if the vector size is less than 2.
 *
 */

Double_t y_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p){

    if (p.size()>= 2){return abs((p[0] + p[1]).Rapidity());}

    return 0;

}

/**
 * @brief Computes the optimized angle of the sum of the first two four-momenta
 *        in the input vector.
 *
 * @param p Vector of four-momenta (e.g. TLorentzVector) representing
 *          reconstructed or generated muons.
 *
 * @return Optimized angle of p[0] + p[1]. Returns 0.0 if the vector size is less than 2.
 *
 */

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

/**
 * @brief Determines whether an event contains a generated Z⁰ decay into two muons, based on particle PDG IDs and mother–daughter relationships.
 *
 * @param n_g Number of generated particles.
 * @param id Vector of PDG IDs for generated particles.
 * @param mid Vector of indices of the mother particle for each generated particle.
 *            The index refers to the position in the `id` vector.
 *
 * @return true if the event contains a Z⁰ (PDG ID = 23) decaying into the
 *         target final state (e.g. muon pair), with consistent mother–daughter
 *         matching; false otherwise.
 */

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

/**
 * @brief Determines whether an event is compatible with a reconstructed Z⁰ decay into a muon pair based on kinematic and isolation criteria.
 *
 * @param n Number of reconstructed particles.
 * @param charge Vector of particle charges.
 * @param Iso Vector of relative isolation values for each particle.
 * @param pt Vector of transverse momentum (pT).
 * @param eta Vector of pseudorapidity (η).
 *
 * @return true if the event contains at least two leptons forming a valid Z⁰ candidate,
 *         typically requiring:
 *         - opposite charge pair,
 *         - transverse momentum above 25 GeV,
 *         - isolation below 0.15,
 *         - pseudorapidity within 2.4 and -2.4;
 *         false otherwise.
 */

bool IsReco(const UInt_t &n, const ROOT::RVec<int> &charge, const ROOT::RVec<float> &Iso,const ROOT::RVec<float> &pt, const ROOT::RVec<float> &eta){
    
    bool IsReco = false;

    if (n == 2){
        if (charge[0] + charge[1] == 0 && Iso[0] < 0.15 && Iso[1] < 0.15 && pt[0] > 25 && pt[1] > 25 && abs(eta[0]) < 2.4 && abs(eta[1]) < 2.4){

            IsReco = true;

        }
    }

    return IsReco;
}

/**
 * @brief Computes four-momenta for selected particles in a generated event,
 *        combining generator-level information with reconstructed kinematics.
 *
 * @param n Number of generated particles.
 * @param id Vector of PDG IDs for generated particles.
 * @param mid Vector of indices of the mother particle (same indexing as `id`).
 * @param pt Vector of transverse momentum (pT) for reconstructed particles.
 * @param eta Vector of pseudorapidity (η).
 * @param phi Vector of azimuthal angle (φ).
 * @param mass Vector of particle masses.
 *
 * @return Vector of four-momenta of particles that pass the selection criteria.
 */

ROOT::RVec<ROOT::Math::PtEtaPhiMVector> GenSel(const UInt_t &n, const ROOT::RVec<int> &id, const ROOT::RVec<int> &mid, const ROOT::RVec<float> &pt, const ROOT::RVec<float> &eta, const ROOT::RVec<float> &phi, const ROOT::RVec<float> &mass){

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

/**
 * @brief Builds four-momenta for selected reconstructed particles in events compatible with a generated Z⁰ decay.

 * @param n_reco Number of reconstructed particles.
 * @param charge Vector of particle charges.
 * @param Iso Vector of relative isolation values.
 * @param pt Vector of transverse momentum (pT).
 * @param eta Vector of pseudorapidity (η).
 * @param phi Vector of azimuthal angle (φ).
 * @param mass Vector of particle masses.
 *
 * @return Vector of four-momenta of reconstructed particles passing the selection criteria.
 */

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

/**
 * @brief Checks whether an invariant mass value falls within a predefined range.
 *
 * @param m Invariant mass (typically in GeV).
 *
 * @return true if m is within the selected mass window, false otherwise.
 *
 * @note The mass window is defined internally in the function, is |m - 91.1817| < 15.
 */

bool Minv_Range(Double_t m) {
    if (std::abs(m - 91.1817) < 15) {return true;}
    return false;
}
