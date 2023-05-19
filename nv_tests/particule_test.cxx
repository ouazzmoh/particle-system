// Include necessary headers
#include <gtest/gtest.h>
#include <particule.hxx>

// Common setup code for all tests
class SetUpForAllTests : public ::testing::Test {
protected:
    Particle p1;
    Particle p2;
    void SetUp() override {
        // Common setup code for all tests
        Vecteur position = Vecteur(3,2,1);
        Vecteur vitesse = Vecteur(1,1,1);
        double masse = 1;
        int identifiant = 2;
        Vecteur force = Vecteur(1,1,1);
        string type = "(nil)";
        p1 = Particle(position, vitesse, masse, identifiant, force, type);
        Vecteur position_2 = Vecteur(5,2.5,4);
        p2 = Particle(position_2, vitesse, masse, identifiant, force, type);
    }

    void TearDown() override {
        // Common teardown code for all tests
    }
};


// On test que l'opérateur < 
TEST_F(SetUpForAllTests, HandlesInfOperator) {
    EXPECT_TRUE(p1 < p2);
}

// On test que l'opérateur ==
TEST_F(SetUpForAllTests, HandlesInfOperator) {
    EXPECT_TRUE(p1 ==  p1);
}
