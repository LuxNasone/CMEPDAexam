#ifndef HISTSTRUCT_H
#define HISTSTRUCT_H

struct OutH{

    std::vector<TH1D*> out;

    void Fill(std::initializer_list<TH1D*> hs);
    
};

#endif