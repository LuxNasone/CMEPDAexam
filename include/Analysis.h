#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <TH1D.h>
#include <vector>
#include <string>

/**
*@brief Block that reconstructs the distributions of variables used to express the differential cross-section (Z transverse momentum, rapidity and optimized angle) from data, 
*       without any unfolding applied. The function takes input data from a folder containing .root files in NANOAD format, processes the events, and returns three histograms (TH1D) that are used for further analysis.
*@param folder_name Path to the input folder cointaining ROOT file with Montecarlo data. File are expected in the NANOAD format.
*@param outname Name for the output file, advised to use a .root file in order to easily see results with a TBrowser. 
*               Default : "NotUnfolded.root";
*@param MT Bool, if true enables multithreading with the option ROOTEnableImplicitMT(). 
*          Default : true;
*@param mute Bool that if true disables some secondary comments. 
*            Default : false.
*@return Vector of the three histograms (TH1D) saved on the file called outname, meant to be unfolded later on.
*        In order : 
*         - [0] : Z tranverse momentum;
*         - [1] : Optimzed angle;
*         - [2] : Z rapidity absolute value;
*        File are saved and can be visualized with a TBrowser.
*@note Requires ROOT framework. 
*@warning Input file must have expected NanoAD format for CMS OpenData.
*/

std::vector<TH1D> NotUnfolded(const char* folder_name, const char* outname = "Repo/outFiles/NotUnfolded.root", bool MT = true, bool mute = false);

/**
* @brief Block that calculates the  response matrices by matching generated and reconstructed events of interest (Z in dimuon).
*        The function builds three RooUnfoldResponse objects (one per observable) and stores them in an output ROOT file, to be used for further analysis.
* @param folder_name Path to the input folder cointaining ROOT file with Montecarlo data. Files are expected in the NANOAD format.
* @param outname Path to the output ROOT file where the response matrices are saved.
*                Default: "Response.root".
* @param MT If true, enables ROOT implicit multithreading via ROOT::EnableImplicitMT().
*           Default: true.
* @param mute Bool that if true disables some secondary comments.
*             Default: false.
*
* @return Void. The function writes RooUnfoldResponse objects and related histograms to the output ROOT file, accessible via TBrowser.
*
* @note Requires ROOT framework and RooUnfold package.
*
* @warning Input file must have expected NanoAD format for CMS OpenData.
*/

void Response(const char* fname, const char* outname = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Response.root", bool MT = true, bool mute = false);

/**
* @brief CrossSection wrapper that applies Bayesian unfolding (RooUnfoldBayes) to three reconstructed histograms, using a precomputed response matrix.
*
* @param folder_name : Path to the input folder cointaining ROOT file with Montecarlo data. Files are expected in the NANOAD format;
* @param n_iter Number of iterations for the Bayesian unfolding algorithm.
* @param rpath Path to the ROOT file containing the response matrix (generated with Response.cpp). The matrices must be compatible in binning and observable definition.
               Default : "Response.root"
* @param outname Name of the output ROOT file where unfolded histograms will be stored.
*                Default : "Unfolded.root"
*
* @return Void. The function writes Unfolded histograms to the output ROOT file, accessible via TBrowser.
*
* @note Requires ROOT framework and RooUnfold package.
*
* @warning Input file must have expected NanoAD format for CMS OpenData.
*/

void Unfolded(const char* folder_name, int n_iter, const char* rpath = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Response.root", const char* outname = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Unfolded.root");

/**
* @brief Block that compares not unfolded and unfolded distributions by overlaying the corresponding histograms. 
*        Also makes a graphs of normalized measured and simulated distributions.
*        Mainly intended to make stylish graphs and for diagnostic of previous outputs.
*
* @param f1 : path to file .root containing not unfolded histograms, assumed to be output of NotUnfolded;
* @param f2 : path to file .root containing unfolded histograms, assumed to be output of Unfolded;
* @param f3 : path to file .root contaning MC histograms, assumed to be output of Response;
* @param outname Name of the output ROOT file where comparison plots are saved.
*                Default : "Comparison.root".
*
* @return Void. The function produces comparison plots (TCanvas objects) stored in the output file and accessible via TBrowser.
*
* @note Requires ROOT framework.
*
* @warning Input file must be output of respectively CrossSection, Response and Unfolded.
*/

void Comparison(const char* f1, const char* f2, const char* f3, const char* outname = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Comparison.root", bool mute = false);

#endif