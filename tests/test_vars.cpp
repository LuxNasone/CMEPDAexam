#include <gtest/gtest.h>

#include "Vars.h"
#include <ROOT/RVec.hxx>
#include <RtypesCore.h>  
#include <Math/Vector4D.h>

TEST(MinvTest, TwoMuons) {

    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> v;

    v.push_back({50, 0.1, 0.2, 0.105});
    v.push_back({40, -0.1, -0.2, 0.105});

    double m = Minv_calculator(v);

    EXPECT_GT(m, 0);
}

TEST(PtTest, SumPt) {

    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> v;

    v.push_back({30, 0, 0, 0.105});
    v.push_back({20, 0, 0, 0.105});

    EXPECT_DOUBLE_EQ(Pt_calculator(v), 50);
}

TEST(RapidityTest, Basic) {
    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> v;
    v.push_back({50, 1.0, 0.0, 0.105});
    v.push_back({50, -1.0, 0.0, 0.105});

    EXPECT_GE(y_calculator(v), 0.0); 
}

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

TEST(MassWindow, Inside) {
    EXPECT_TRUE(Minv_Range(91.2));
    EXPECT_FALSE(Minv_Range(120.0));
}