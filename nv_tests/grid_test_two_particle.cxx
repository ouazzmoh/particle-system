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
    Universe *univ;
    Particle* particleToInsert;
    Particle* secondParticle;
    void SetUp() override {
        // Common setup code for all tests
        double masse = 1;
        string type = "(nil)";
        Vecteur force = Vecteur(0.0, 0.0, 0.0);
        Vecteur position = Vecteur(0, 0,0);
        Vecteur vitesse =  Vecteur(7.0, 7.0, 7.0);
        particleToInsert = new Particle(position, vitesse, masse, 1, force, type);

        vector<Particle*> particleVec;
        particleVec.push_back(particleToInsert);

        int dim = 2;
        double epsilon = 1;
        double rCut = 10;
        double sigma = 1;
        double LD[2] = {20,20};
        int boundCond =0;

        // we create a second particle
        double masse_2 = 1;
        string type_2 = "(nil)";
        Vecteur force_2 = Vecteur(0.0, 0.0, 0.0);
        Vecteur position_2 = Vecteur(15, 4 ,9);
        Vecteur vitesse_2 =  Vecteur(-3.0, 0.0, 5.0);
        secondParticle = new Particle(position_2, vitesse_2, masse_2, 2, force_2, type_2);

        particleVec.push_back(secondParticle);
        univ = new Universe(particleVec, dim, rCut, LD, epsilon, sigma, boundCond);
    }

    void TearDown() override {
        // Common teardown code for all tests
        delete univ;
        vector<Particle*> p = univ->getParticles();
        for (Particle* particle : p) {
            delete particle;
        }
        delete particleToInsert;
        delete secondParticle;
    }
};







// check the universe after the application of stromerVerletPotential function
// when there is 2 particles
TEST_F(SetUpForAllTests, Handles2Particles) {
    // The expected position for the second particle is grid[1][0]
    unordered_set<Particle*> particlesIndex_zero_zero = univ->getGrid()[0][0].getParticles();
    bool b  = checkParticleInCell(secondParticle, particlesIndex_zero_zero);
    EXPECT_TRUE(!b);

    unordered_set<Particle*> particlesIndex_zero_one = univ->getGrid()[0][1].getParticles();
    bool c  = checkParticleInCell(secondParticle, particlesIndex_zero_one);
    EXPECT_TRUE(!c);

    unordered_set<Particle*> particlesIndex_one_zero = univ->getGrid()[1][0].getParticles();
    bool d  = checkParticleInCell(secondParticle, particlesIndex_one_zero);
    EXPECT_TRUE(d);

    unordered_set<Particle*> particlesIndex_one_one = univ->getGrid()[1][1].getParticles();
    bool e  = checkParticleInCell(secondParticle, particlesIndex_one_one);
    EXPECT_TRUE(!e);

    
    // we applicate stromerVerletPotential
    stromerVerletPotential(*univ, 2,0.5, false, "", false, 0, 0);


    // The new expected position for second particle is grid[0][0]
    // because it has a negative velocity on the direction x and a 0 velocity in the direction y
    // The new expected position for first particle is grid[1][1]
    // because it has a positive velocity on the direction x and y

    // We check that for the first particle
    particlesIndex_zero_zero = univ->getGrid()[0][0].getParticles();
    b  = checkParticleInCell(particleToInsert, particlesIndex_zero_zero);
    EXPECT_TRUE(!b);

    particlesIndex_zero_one = univ->getGrid()[0][1].getParticles();
    c  = checkParticleInCell(particleToInsert, particlesIndex_zero_one);
    EXPECT_TRUE(!c);

    particlesIndex_one_zero = univ->getGrid()[1][0].getParticles();
    d  = checkParticleInCell(particleToInsert, particlesIndex_one_zero);
    EXPECT_TRUE(!d);

    particlesIndex_one_one = univ->getGrid()[1][1].getParticles();
    e  = checkParticleInCell(particleToInsert, particlesIndex_one_one);
    EXPECT_TRUE(e);


    // We check that for the second particle
    particlesIndex_zero_zero = univ->getGrid()[0][0].getParticles();
    b  = checkParticleInCell(secondParticle, particlesIndex_zero_zero);
    EXPECT_TRUE(b);

    particlesIndex_zero_one = univ->getGrid()[0][1].getParticles();
    c  = checkParticleInCell(secondParticle, particlesIndex_zero_one);
    EXPECT_TRUE(!c);

    particlesIndex_one_zero = univ->getGrid()[1][0].getParticles();
    d  = checkParticleInCell(secondParticle, particlesIndex_one_zero);
    EXPECT_TRUE(!d);

    particlesIndex_one_one = univ->getGrid()[1][1].getParticles();
    e  = checkParticleInCell(secondParticle, particlesIndex_one_one);
    EXPECT_TRUE(!e);
}


