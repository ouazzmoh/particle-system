// Include necessary headers
#include <gtest/gtest.h>
#include "vecteur.hxx"
#include "universe.hxx"



using namespace std;

// Common setup code for all tests
class SetUpForAllTests : public ::testing::Test {
    protected:
        Vecteur v1;
        Vecteur v2;
        void SetUp() override {
            // Common setup code for all tests
            v1 = Vecteur(3.0, 2.0, 1.0);
            v2 = Vecteur(2,1,1);
        }

        void TearDown() override {
            // Common teardown code for all tests
        }
};


// On test que la norme de 0 est 0
TEST_F(SetUpForAllTests, HandlesZeroInput) {
    Vecteur v = Vecteur(0, 0, 0);
    double expectedNorm = 0;
    EXPECT_EQ(v.norm(), 0);
}


// On teste que la norme est positive 
TEST_F(SetUpForAllTests, HandlesNormPositivity) {
    v2 = Vecteur(1, -4, 8);
    Vecteur v3 = Vecteur(-7, -8, -4);
    EXPECT_GT(v1.norm(), 0);
    EXPECT_GT(v2.norm(), 0);
    EXPECT_GT(v3.norm(), 0);
}

// On teste que la norme d'un vecteur est très proche de la vrai valeur
// de la norme caculée séparément
TEST(VecteurNormTest, HandleNormOp) {
    Vecteur v3 = Vecteur(1, 1, 1);
    double expectedNorm = 1.73205080757;
    double tolerance = 0.01; // Define the tolerance for comparison
    EXPECT_NEAR(v3.norm(), expectedNorm, tolerance);
}


// On teste l'opérateur + 
TEST_F(SetUpForAllTests, HandlesPlusOperator) {
    Vecteur v3 = Vecteur(-2,-1,-1);
    EXPECT_TRUE(v1+v2 ==  Vecteur(5,3,2));
    EXPECT_TRUE(v3+v2 ==  Vecteur(0,0,0));
}

// On teste l'opérateur += 
TEST_F(SetUpForAllTests, HandlesPlusEqualOperator) {
    v1 += v2;
    EXPECT_TRUE(v1 ==  Vecteur(5,3,2));
}


// On teste l'opérateur *
TEST_F(SetUpForAllTests, HandlesMultOperator) {
    EXPECT_TRUE(v1*v2 ==  Vecteur(6,2,1));
}

// On teste l'opérateur *=
TEST_F(SetUpForAllTests, HandlesMultEqOperator) {
    v1 *= v2;
    EXPECT_TRUE(v1 ==  Vecteur(6,2,1));
}

// On teste l'opérateur -
TEST_F(SetUpForAllTests, HandlesMinusOperator) {
    Vecteur v1 = Vecteur(3.0, 2.0, 1.0);
    Vecteur v2 = Vecteur(2,1,1);
    EXPECT_TRUE(v1-v2 ==  Vecteur(1,1,0));
}

// On teste l'opérateur -=
TEST_F(SetUpForAllTests, HandlesMinusEqOperator) {
    Vecteur v1 = Vecteur(3.0, 2.0, 1.0);
    Vecteur v2 = Vecteur(2,1,1);
    v1 -= v2;
    EXPECT_TRUE(v1 ==  Vecteur(1,1,0));
}

// On teste l'opérateur ==
TEST_F(SetUpForAllTests, HandlesEqEqOperator) {
    v2 = Vecteur(3, 2, 1);
    EXPECT_TRUE(v1 == v2);
}


