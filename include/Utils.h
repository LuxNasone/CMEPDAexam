#ifndef UTILS_H
#define UTILS_H

#include <Rtypes.h>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>

/** @defgroup GlobalVariables Global Variables */

/// @ingroup GlobalVariables
/// @brief Number of bins, used for all histograms and the response matrices. 

extern int n_b;

/// @ingroup GlobalVariables
/// @brief Variable names, used in loops to define new RDataFrame columns 

extern std::vector<const char*> vars;

/// @ingroup GlobalVariables
/// @brief X axis names, used in loops for plots, only meant for aesthetic usage.

extern std::vector<const char*> xlabels;

/// @ingroup GlobalVariables
/// @brief Y axis names, used in loops for plots, only meant for aesthetic usage in the final cross section plot.

extern std::vector<const char*> ylabels;

/// @ingroup GlobalVariables
/// @brief Title for graphs, used in loops for plots, only meant for aesthetic usage.

extern std::vector<const char*> titles;

/// @ingroup GlobalVariables
/// @brief Bounds for variables, both for graphs and ranges in response matrices. 

extern std::vector<std::pair<Float_t, Float_t>> bounds;

/// @ingroup GlobalVariables
/// @brief Bounds for y axis, purely aesthetic in unfolded graphs.

extern std::vector<std::pair<Float_t, Float_t>> range;

/// @ingroup GlobalVariables
/// @brief Integrated luminosity for used dataset (expressed in [1/pb]).

extern double L;

/**
* @brief Computes the invariant mass of the sum of the first two four-momentain the input vector.
*
* @param p Vector of four-momenta (PtEtaPhiMVector) representing reconstructed or generated muons.
*
* @return Invariant mass of p[0] + p[1]. Returns 0.0 if the vector size is less than 2.
*
*/

Double_t Minv_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p);

/**
* @brief Computes the transverse mmomentum of the sum of the first two four-momenta in the input vector.
*
* @param p Vector of four-momenta (PtEtaPhiMVector) representing reconstructed or generated muons.
*
* @return Transverse momentum of p[0] + p[1]. Returns 0.0 if the vector size is less than 2.
*
*/

Double_t Pt_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p);

/**
* @brief Computes the rapidity absolute value of the sum of the first two four-momenta in the input vector.
*
* @param p Vector of four-momenta (PtEtaPhiMVector) representing reconstructed or generated muons.
*
* @return Rapidity absolute value of p[0] + p[1]. Returns 0.0 if the vector size is less than 2.
*
*/

Double_t y_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p);

/**
* @brief Computes the optimized angle of the two four-momenta in the input vector.
*
* @param p Vector of four-momenta (PtEtaPhiMVector) representing reconstructed or generated muons.
*
* @return Optimized angle of p[0] and p[1]. Returns 0.0 if the vector size is less than 2.
*
*/

Double_t phi_eta_calculator(const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p);

/**
* @brief Computes four-momenta for selected particles in a generated event, combining generator-level information with reconstructed kinematics.
*
* @param n Number of generated particles.
* @param id Vector of PDG IDs for generated particles.
* @param mid Vector of indices of the mother particle for each generated particle. The index refers to the position in the id vector.
* @param pt Vector of transverse momentums.
* @param eta Vector of pseudorapidities.
* @param phi Vector of azimuthal angles.
* @param mass Vector of particle masses.
*
* @return Vector of four-momenta of particles that pass the selection criteria.
*/

ROOT::RVec<ROOT::Math::PtEtaPhiMVector> GenSel(const UInt_t &n, const ROOT::RVec<int> &id, const ROOT::RVec<int> &mid, const ROOT::RVec<float> &pt, const ROOT::RVec<float> &eta, const ROOT::RVec<float> &phi, const ROOT::RVec<float> &mass);

/**
* @brief Builds four-momenta for selected reconstructed particles in events compatible with a generated Z⁰ decay.
*
* @param n_reco Number of reconstructed particles.
* @param charge Vector of particle charges.
* @param Iso Vector of relative isolation values.
* @param pt Vector of transverse momentums .
* @param eta Vector of pseudorapidities.
* @param phi Vector of azimuthal angles.
* @param mass Vector of particle masses.
*
* @return Vector of four-momenta of reconstructed particles passing the selection criteria.
*/

ROOT::RVec<ROOT::Math::PtEtaPhiMVector> Reco(const UInt_t &n, const ROOT::RVec<int> &charge, const ROOT::RVec<float> &Iso,const ROOT::RVec<float> &pt, const ROOT::RVec<float> &eta, const ROOT::RVec<float> &phi,const ROOT::RVec<float> &mass);

/**
* @brief Checks whether an invariant mass value falls within a predefined range.
*
* @param m Invariant mass (expressed in GeV).
*
* @return true if m is within the selected mass window, false otherwise.
*
* @note The mass window is defined internally in the function, is |m - 91.1817| < 15.
*/

bool Minv_Range(Double_t m);

/**
* @brief Determines whether an event is compatible with a reconstructed Z decay into a muon pair based on kinematic and isolation criteria.
*
* @param n Number of reconstructed particles.
* @param charge Vector of particle charges.
* @param Iso Vector of relative isolation values for each particle.
* @param pt Vector of transverse momentums.
* @param eta Vector of pseudorapidities.
* @param mass Vector of particle masses.
*
* @return true if the event contains at least two leptons forming a valid Z candidate, requiring:
*         - opposite charge pair;
*         - transverse momentum above 25 GeV;
*         - isolation below 0.15;
*         - pseudorapidity within 2.4 and -2.4;
*         - |mass -0.1057| < 2.5e-5
*         false otherwise.
*/

bool IsReco(const UInt_t &n, const ROOT::RVec<int> &charge, const ROOT::RVec<float> &Iso,const ROOT::RVec<float> &pt, const ROOT::RVec<float> &eta, const ROOT::RVec<float> &mass);

/**
* @brief Determines whether an event contains a generated Z decay into two muons, based on particle PDG IDs and mother–daughter relationships.
*
* @param n_g Number of generated particles.
* @param id Vector of PDG IDs for generated particles.
* @param mid Vector of indices of the mother particle for each generated particle. The index refers to the position in the id vector.
*
* @return true if the event contains a Z (PDG ID = 23) decaying into the target final state (e.g. muon pair), with consistent mother–daughter matching; false otherwise.
*/

bool IsTrue(const UInt_t &n_g,const ROOT::RVec<int> &id, const ROOT::RVec<int> &mid);

#endif