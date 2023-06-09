#include "universe.hxx"
#include "particule.hxx"
#include <random>
#include <gtest/gtest.h>



using namespace std;

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
    Universe *univ;
    Particle* particleToInsert;
    void SetUp() override {
        // Common setup code for all tests
        double masse = 1;
        string type = "(nil)";
        Vecteur force = Vecteur(0.0, 0.0, 0.0);
        Vecteur position = Vecteur(10, 10,0);
        Vecteur vitesse =  Vecteur(0.0, 0.0, 7.0);
        particleToInsert = new Particle(position, vitesse, masse, 1, force, type);

        vector<Particle*> particleVec;
        particleVec.push_back(particleToInsert);

        int dim = 2;
        double epsilon = 1;
        double rCut = 10;
        double sigma = 1;
        double LD[2] = {20,20};
        int boundCond =0;
        univ = new Universe(particleVec, dim, rCut, LD, epsilon, sigma, -1);
    }

    void TearDown() override {
        // Common teardown code for all tests
        delete univ;
        vector<Particle*> p = univ->getParticles();
        for (Particle* particle : p) {
            delete particle;
        }
        delete particleToInsert;
    }
};

// check the universe after the application of stromerVerletPotential function
// in presence of a gravitational field
TEST_F(SetUpForAllTests, HandlesStormerVerlet) {
    stromerVerletPotential(*univ, 2,0.5, false, "", false, 9.8,0);
    // the initial position of the particle is grid[1][1]
    // The new expected position for the particle is grid[1][0]
    // because of the gravitational field ( the particle is supposed to fall )
    unordered_set<Particle*> particlesIndex_zero_zero = univ->getGrid()[0][0].getParticles();
    bool b  = checkParticleInCell(particleToInsert, particlesIndex_zero_zero);
    EXPECT_TRUE(!b);

    unordered_set<Particle*> particlesIndex_zero_one = univ->getGrid()[0][1].getParticles();
    bool c  = checkParticleInCell(particleToInsert, particlesIndex_zero_one);
    EXPECT_TRUE(!c);

    unordered_set<Particle*> particlesIndex_one_zero = univ->getGrid()[1][0].getParticles();
    bool d  = checkParticleInCell(particleToInsert, particlesIndex_one_zero);
    EXPECT_TRUE(d);

    unordered_set<Particle*> particlesIndex_one_one = univ->getGrid()[1][1].getParticles();
    bool e  = checkParticleInCell(particleToInsert, particlesIndex_one_one);
    EXPECT_TRUE(!e);
    
}



