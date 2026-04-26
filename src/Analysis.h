#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <TH1D.h>
#include <vector>
#include <string>

std::vector<TH1D> CrossSection(const char* folder_name, std::string outname = "Repo/outFiles/NotUnfolded.root", bool MT = true, bool mute = false);

void Response(const char* fname, std::string outname = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Response.root", bool MT = true, bool mute = false);

void Unfolded(const char* folder_name, int n_iter, const char* rpath = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Response.root", const char* outname = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Unfolded.root");

void Comp(const char* f1, const char* f2, const char* f3, const char* outname = "/home/lux_n/CMEPDA/Exam/Repo/outFiles/Comparison.root", bool mute = false);

#endif