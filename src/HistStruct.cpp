#include <TH1D.h>
#include <vector> 
#include "HistStruct.h"

void OutH::Fill(std::initializer_list<TH1D*> hs){

    for(auto h : hs){out.push_back(h);}

};
