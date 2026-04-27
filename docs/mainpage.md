## Differential Z boson production cross section in proton-proton collision at 13 TeV with CMS Open Data

Repository on the project made for the exam of Computing Methods for Experimental Physics and Data Analysis.

## Description

The task is to reproduce the measurement of the differential production cross section for Z in p-p collision at 13 TeV with Open Data fournished by CMS.

For the article we refer to: https://cms-results.web.cern.ch/cms-results/public-results/publications/SMP-17-010/index.html.

We select muon pairs applying the cuts described in the article:

- n = 2;
- pt_{i}> 25 GeV i=1,2;
- |η_{i}| < 2.4 i=1,2;
- ΔR_{i} < 0.15 i=1,2;
- q_1 +q_2 = 0;

Finding the usual Breight-Wigner peak around expected Z mass. After finding the resonance of the Z and selecting events with |m_{inv}−m{Z^{0}}|< 15 GeV, we calculate distributions for the following quantities:

- Transverse momentum : P_{t};
- Rapidity : y;
- Optimized angle : ϕ_{η} = tan((π−Δϕ)/2)sin(θ_{η});

With cos(θ_{η})=tanh(Δη/2).

In order to to take into account inefficiencies we performed an unfolding procedure, estimating the response matrix with data generated from a Montecarlo and using the RooUnfold toolbox. The generated events are determined to be Z decays in two muons by:

- Selecting ±13 PDG ID to select muon pairs;
- Finds index for mother particle and controls if it is 23 (Z PDG ID);
Then on the same dataset we apply the cut of the article to have a generated-reconstructed match. The response matrix is then obtained by using RooUnfoldResponse. We then use a RooUnfoldBayesian oon measured distribution.

The differential cross section is obtain by scaling the histograms with the integrated luminosity and bin width. The dataset used has an integrated luminosity 
L = 8740 pb^-1.

## How to use

To compile the program one must : 

- Use ROOT (type root or root -l on linux shell);
- Import the RooUnfold package with gSystem->Load("/path/libRooUnfold.so");
- Compile Vars.cpp (.L /path/Vars.cpp+);
- Compile Analysis.cpp (.L /path/Analysis.cpp+);

Then you can execute all code blocks in there. The typical workflow is to run them sequentially to generate all .root file and then compare them with Comp script, which requires all three output.
