#include <gtest/gtest.h>

#include <Main.h>
#include <ROOT/RVec.hxx>
#include <RtypesCore.h>  
#include <Math/Vector4D.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <vector>

//Global variables testing

extern int n_b;
extern std::vector<std::string> vars;
extern std::vector<std::pair<Float_t, Float_t>> bounds;
extern std::vector<std::string> xlabels;
extern std::vector<std::string> ylabels;

//Mismatch test

TEST(GlobalConfigTest, SizesAreConsistent) {
    EXPECT_EQ(vars.size(), bounds.size());
    EXPECT_EQ(vars.size(), xlabels.size());
    EXPECT_EQ(vars.size(), ylabels.size());
}

//Bounds min < max consistency

TEST(GlobalConfigTest, BoundsAreValid) {
    for (const auto& b : bounds) { EXPECT_LT(b.first, b.second);}
}

//This block below creates a test .root file

void TestFile(const char* filename = "test_data.root") {

    TFile file(filename, "RECREATE");
    TTree tree("Events", "Test Tree");

    int nMuon;
    std::vector<float> Muon_pt, Muon_eta, Muon_phi, Muon_mass, Muon_pfRelIso03_all;
    std::vector<int> Muon_charge;

    tree.Branch("nMuon", &nMuon);
    tree.Branch("Muon_pt", &Muon_pt);
    tree.Branch("Muon_eta", &Muon_eta);
    tree.Branch("Muon_phi", &Muon_phi);
    tree.Branch("Muon_mass", &Muon_mass);
    tree.Branch("Muon_pfRelIso03_all", &Muon_pfRelIso03_all);
    tree.Branch("Muon_charge", &Muon_charge);

    nMuon = 2;
    Muon_pt = {30, 28};
    Muon_eta = {0.5, -0.3};
    Muon_phi = {0.1, 3.0};
    Muon_mass = {0.105, 0.105};
    Muon_pfRelIso03_all = {0.1, 0.1};
    Muon_charge = {+1, -1};
    tree.Fill();

    nMuon = 2;
    Muon_pt = {35, 40};
    Muon_eta = {0.2, -0.2};
    Muon_phi = {1.0, 2.0};
    Muon_mass = {0.105, 0.105};
    Muon_pfRelIso03_all = {0.2, 0.1};
    Muon_charge = {+1, -1};
    tree.Fill();

    nMuon = 2;
    Muon_pt = {50, 45};
    Muon_eta = {0.1, -0.1};
    Muon_phi = {0.5, 2.5};
    Muon_mass = {0.105, 0.105};
    Muon_pfRelIso03_all = {0.1, 0.1};
    Muon_charge = {+1, +1}; // FAIL
    tree.Fill();

    nMuon = 2;
    Muon_pt = {20, 30}; // FAIL
    Muon_eta = {0.3, -0.3};
    Muon_phi = {0.2, 3.1};
    Muon_mass = {0.105, 0.105};
    Muon_pfRelIso03_all = {0.1, 0.1};
    Muon_charge = {+1, -1};
    tree.Fill();

    nMuon = 2;
    Muon_pt = {30, 30};
    Muon_eta = {3.0, 0.1}; 
    Muon_phi = {0.2, 3.1};
    Muon_mass = {0.105, 0.105};
    Muon_pfRelIso03_all = {0.1, 0.1};
    Muon_charge = {+1, -1};
    tree.Fill();

    file.Write();
    file.Close();
}

//Cross section test : correct output shape

TEST(CrossSectionTest, ReturnsThreeHistograms) {

    auto result = CrossSection("test_data", "test.root", false, true);

    EXPECT_EQ(result.size(), 3);

}

//Cross section test : correct number of bins

TEST(CrossSectionTest, HistogramBinsAreCorrect) {

    auto result = CrossSection("test_data", "test.root", false, true);

    for (const auto& h : result) { EXPECT_EQ(h.GetNbinsX(), n_b); }

}

