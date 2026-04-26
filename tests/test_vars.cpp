#include <gtest/gtest.h>

#include "Vars.h"
#include <ROOT/RVec.hxx>
#include <RtypesCore.h>  
#include <Math/Vector4D.h>

//Minv_calculator test

TEST(MinvTest, TwoMuons) {

    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> v;

    v.push_back({50, 0.1, 0.2, 0.105});
    v.push_back({40, -0.1, -0.2, 0.105});

    double m = Minv_calculator(v);

    EXPECT_GT(m, 0);
}

TEST(MinvTest, LessThanTwo) {
    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> v;
    v.push_back({50, 0, 0, 0.105});

    EXPECT_DOUBLE_EQ(Minv_calculator(v), 0.0);
}

//Pt_calculator test

TEST(PtTest, SumPt) {

    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> v;

    v.push_back({30, 0, 0, 0.105});
    v.push_back({20, 0, 0, 0.105});

    EXPECT_DOUBLE_EQ(Pt_calculator(v), 50);
}

TEST(PtTest, LessThanTwo) {
    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> v;

    EXPECT_DOUBLE_EQ(Pt_calculator(v), 0.0);
}

//y_calculator test

TEST(RapidityTest, Basic) {
    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> v;
    v.push_back({50, 1.0, 0.0, 0.105});
    v.push_back({50, -1.0, 0.0, 0.105});

    EXPECT_GE(y_calculator(v), 0.0); 
}

TEST(RapidityTest, LessThanTwo) {
    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> v;

    EXPECT_DOUBLE_EQ(y_calculator(v), 0.0);
}

//phi_eta_calculator tests

TEST(PhiEtaTest, Basic) {
    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> v;
    v.push_back({50, 0.1, 0.0, 0.105});
    v.push_back({50, -0.1, 0.0, 0.105});

    EXPECT_TRUE(std::isfinite(phi_eta_calculator(v)));
}

TEST(PhiEtaTest, LessThanTwo) {
    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> v;
    v.push_back({50, 0.1, 0.1, 0.105});

    EXPECT_DOUBLE_EQ(phi_eta_calculator(v), 0.0);
}

//IsTrue test

TEST(IsTrueTest, ValidZToMuMu) {
    ROOT::RVec<int> id  = {23, 13, -13};
    ROOT::RVec<int> mid = {-1, 0, 0};

    EXPECT_TRUE(IsTrue(3, id, mid));
}

TEST(IsTrueTest, SameChargeFail) {
    ROOT::RVec<int> id  = {23, 13, 13};
    ROOT::RVec<int> mid = {-1, 0, 0};

    EXPECT_FALSE(IsTrue(3, id, mid));
}

TEST(IsTrueTest, WrongMotherFail) {
    ROOT::RVec<int> id  = {24, 13, -13};
    ROOT::RVec<int> mid = {-1, 0, 0};

    EXPECT_FALSE(IsTrue(3, id, mid));
}

TEST(IsTrueTest, LessThanTwoMuons) {
    ROOT::RVec<int> id  = {23, 13};
    ROOT::RVec<int> mid = {-1, 0};

    EXPECT_FALSE(IsTrue(2, id, mid));
}

//GenSel test

TEST(GenSelTest, ValidZToMuMu) {

    ROOT::RVec<int> id  = {23, 13, -13};
    ROOT::RVec<int> mid = {-1, 0, 0};

    ROOT::RVec<float> pt   = {0, 40, 45};
    ROOT::RVec<float> eta  = {0, 0.1, -0.1};
    ROOT::RVec<float> phi  = {0, 0.2, -0.2};
    ROOT::RVec<float> mass = {0, 0.105, 0.105};

    auto v = GenSel(3, id, mid, pt, eta, phi, mass);

    EXPECT_EQ(v.size(), 2);
}

TEST(GenSelTest, NoZMotherRejected) {

    ROOT::RVec<int> id  = {24, 13, -13}; // no Z
    ROOT::RVec<int> mid = {-1, 0, 0};

    ROOT::RVec<float> pt   = {0, 40, 45};
    ROOT::RVec<float> eta  = {0, 0.1, -0.1};
    ROOT::RVec<float> phi  = {0, 0.2, -0.2};
    ROOT::RVec<float> mass = {0, 0.105, 0.105};

    auto v = GenSel(3, id, mid, pt, eta, phi, mass);

    EXPECT_EQ(v.size(), 0);
}

TEST(GenSelTest, LessThanTwoMuons) {

    ROOT::RVec<int> id  = {23, 13};
    ROOT::RVec<int> mid = {-1, 0};

    ROOT::RVec<float> pt   = {0, 40};
    ROOT::RVec<float> eta  = {0, 0.1};
    ROOT::RVec<float> phi  = {0, 0.2};
    ROOT::RVec<float> mass = {0, 0.105};

    auto v = GenSel(2, id, mid, pt, eta, phi, mass);

    EXPECT_EQ(v.size(), 0);
}

//IsReco test

TEST(RecoTest, ValidEvent) {

    ROOT::RVec<int> charge = {1, -1};
    ROOT::RVec<float> iso = {0.1, 0.1};
    ROOT::RVec<float> pt = {30, 30};
    ROOT::RVec<float> eta = {1.0, -1.0};

    EXPECT_TRUE(IsReco(2, charge, iso, pt, eta));
}

TEST(RecoTest, IsolationFail) {
    ROOT::RVec<int> charge = {1, -1};
    ROOT::RVec<float> iso = {0.2, 0.1};
    ROOT::RVec<float> pt = {30, 30};
    ROOT::RVec<float> eta = {1.0, -1.0};

    EXPECT_FALSE(IsReco(2, charge, iso, pt, eta));
}

//Reco test

TEST(RecoVecTest, ValidRecoBuildsVectors) {
    ROOT::RVec<int> charge = {1, -1};
    ROOT::RVec<float> iso  = {0.1, 0.1};
    ROOT::RVec<float> pt   = {30, 35};
    ROOT::RVec<float> eta  = {0.5, -0.5};
    ROOT::RVec<float> phi  = {0.1, -0.1};
    ROOT::RVec<float> mass = {0.105, 0.105};

    auto v = Reco(2, charge, iso, pt, eta, phi, mass);

    EXPECT_EQ(v.size(), 2);
}

TEST(RecoVecTest, FailsSelectionReturnsEmpty) {
    ROOT::RVec<int> charge = {1, -1};
    ROOT::RVec<float> iso  = {0.2, 0.1}; 
    ROOT::RVec<float> pt   = {30, 35};
    ROOT::RVec<float> eta  = {0.5, -0.5};
    ROOT::RVec<float> phi  = {0.1, -0.1};
    ROOT::RVec<float> mass = {0.105, 0.105};

    auto v = Reco(2, charge, iso, pt, eta, phi, mass);

    EXPECT_EQ(v.size(), 0);
}

//Minv_Range testy

TEST(MassWindow, Inside) {
    EXPECT_TRUE(Minv_Range(91.2));
    EXPECT_FALSE(Minv_Range(120.0));
}