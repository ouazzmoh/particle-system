#include "universe.hxx"
#include "particule.hxx"
#include <random>
#include <gtest/gtest.h>

bool checkParticleInCell(Particle *p, unordered_set<Particle*> &vecParticles){
    for ( auto &particle : vecParticles){
        if (*particle == *p){
            return true;
        }
    }
    return false;
}

// Common setup code for all tests
class SetUpForAllTests : public ::testing::Test {
    protected:
        Universe *univAbsorption;
        Universe *univReflection;
        Universe *univPeriodic;
        Universe *univRefLJlexion;
        Particle* particleToInsert_1;
        Particle* particleToInsert_2;
        Particle* particleToInsert_3;
        Particle* particleToInsert_4;

        void SetUp() override {
            // Common setup code for all tests
            string type = "nil";

            particleToInsert_1 = new Particle(Vecteur(50.0, 50.0, 0.0),
                                                  Vecteur(5.0, 0.0, 0.0), 1.0, 1, Vecteur(0.0, 0.0, 0.0), type);
       
            vector<Particle*> particleVec_1;
            particleVec_1.push_back(particleToInsert_1);

            particleToInsert_2 = new Particle(Vecteur(50.0, 50.0, 0.0),
                                                  Vecteur(5.0, 0.0, 0.0), 1.0, 1, Vecteur(0.0, 0.0, 0.0), type);
       
            vector<Particle*> particleVec_2;
            particleVec_2.push_back(particleToInsert_2);

            particleToInsert_3 = new Particle(Vecteur(50.0, 50.0, 0.0),
                                                  Vecteur(5.0, 0.0, 0.0), 1.0, 1, Vecteur(0.0, 0.0, 0.0), type);
       
            vector<Particle*> particleVec_3;
            particleVec_3.push_back(particleToInsert_3);

            particleToInsert_4 = new Particle(Vecteur(50.0, 50.0, 0.0),
                                                  Vecteur(5.0, 0.0, 0.0), 1.0, 1, Vecteur(0.0, 0.0, 0.0), type);
       
            vector<Particle*> particleVec_4;
            particleVec_4.push_back(particleToInsert_4);

            int dim = 2;
            double lD_1[2] = {100, 100};
            double lD_2[2] = {100, 100};
            double lD_3[2] = {100, 100};
            double lD_4[2] = {100, 100};

            double epsilon = 5;
            double sigma = 1;
            double rCut = 10;
            univAbsorption = new Universe(particleVec_1, 2, rCut, lD_1, epsilon, sigma, 0);
            univReflection = new Universe(particleVec_2, 2, rCut, lD_2, epsilon, sigma, -1);
            univPeriodic = new Universe(particleVec_3,  2, rCut, lD_3, epsilon, sigma, 1);
            univRefLJlexion = new Universe(particleVec_4, 2, rCut, lD_4, epsilon, sigma, 0);
        }

        void TearDown() override {
            // Common teardown code for all tests
            delete univAbsorption;
            delete univReflection;
            delete univPeriodic;
            delete univRefLJlexion;
        }
};




/*
// test Absorption
TEST_F(SetUpForAllTests, HandlesAbsorption) {
    stromerVerletPotential(*univAbsorption, 20, 0.005, true, "../../demo/abs/", false, 0, 0);
    ADD_FAILURE() << "univ : " << *univAbsorption;
    bool b;
    for( int i = 0; i < univAbsorption->getGrid().size(); i++){
        for( int j = 0; j < univAbsorption->getGrid()[i].size(); j++){
            //ADD_FAILURE() << "i j  : " << i  << "djsh"<< j;
            vector<Particle*> particlesIndex = univReflection->getGrid()[i][j].getParticles();
            b = checkParticleInCell(particleToInsert_1,particlesIndex);
            //ADD_FAILURE() << b;
            if (b  == true){
                ADD_FAILURE() << "i j  : " << i  << "djsh"<< j;
            }
            EXPECT_TRUE(!b);
        } 
    }
}
*/

// test Reflection
TEST_F(SetUpForAllTests, HandlesReflection) {
    // the particle is initially at grid[5][5]
    // and we expect it to return to grid[5][5] after reflection
    stromerVerletPotential(*univReflection, 20, 0.005, true, "../../demo/ref/", false, 0, 0);
    unordered_set<Particle*> particlesIndex_five_five = univReflection->getGrid()[5][5].getParticles();
    bool b  = checkParticleInCell(particleToInsert_2, particlesIndex_five_five);
    EXPECT_TRUE(b);
}


// test periodic
TEST_F(SetUpForAllTests, HandlesPeriodic) {
    // the particle is initially at grid[5][5]
    // and we expect it to return to grid[5][5] when it is periodic
    stromerVerletPotential(*univPeriodic, 20, 0.005, true, "../../demo/per/", false, 0, 0);
    unordered_set<Particle*> particlesIndex_five_five = univReflection->getGrid()[5][5].getParticles();
    bool b  = checkParticleInCell(particleToInsert_2, particlesIndex_five_five);
    EXPECT_TRUE(b);
}


// test RefLJlexion
TEST_F(SetUpForAllTests, HandlesRefLJlexion) {
    // the particle is initially at grid[5][5]
    // and we expect it to return to grid[5][5] when it is periodic
    stromerVerletPotential(*univRefLJlexion, 20, 0.005, true, "../../demo/per/", true, 0, 0);
    double tolerance = 1; // Define the tolerance for comparison
    EXPECT_NEAR(particleToInsert_4->getPosition().getX(),50, tolerance);
    EXPECT_NEAR(particleToInsert_4->getPosition().getY(),50, tolerance);
}


