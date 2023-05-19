#include <universe.hxx>
#include <particule.hxx>
#include <random>
#include <gtest/gtest.h>

bool checkParticleInCell(Particle *p, vector<Particle*> &vecParticles){
    for ( auto &particle : vecParticles){
        if (*particle == *p){
            return true;
        }
    }
    return false;
}

int taillGrid(std::vector<std::vector<Cell>> g){
    int c = 0;
    for (int i = 0; i < g.size(); i++){
        for (int j = 0; j < g[i].size(); j++){
            c++;
        }
    }
    return c;
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
        univ = new Universe(particleVec, dim, rCut, LD, epsilon, sigma, boundCond);
        // On crée un universe qui a donc 4 carrées. 


        // we create a second for the case when we test the universe with 2 particles
        double masse_2 = 1;
        string type_2 = "(nil)";
        Vecteur force_2 = Vecteur(0.0, 0.0, 0.0);
        Vecteur position_2 = Vecteur(15, 4 ,9);
        Vecteur vitesse_2 =  Vecteur(3.0, 0.0, 5.0);
        secondParticle = new Particle(position_2, vitesse_2, masse_2, 1, force_2, type_2);
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

// test the initialization of the grid
TEST_F(SetUpForAllTests, HandlesGridInitialPosition) {

    vector<Particle*> particlesIndex_zero_zero = univ->getGrid()[0][0].getParticles();
    bool b  = checkParticleInCell(particleToInsert, particlesIndex_zero_zero);
    EXPECT_TRUE(b);

    vector<Particle*> particlesIndex_zero_one = univ->getGrid()[0][1].getParticles();
    bool c  = checkParticleInCell(particleToInsert, particlesIndex_zero_one);
    EXPECT_TRUE(!c);

    vector<Particle*> particlesIndex_one_zero = univ->getGrid()[1][0].getParticles();
    bool d  = checkParticleInCell(particleToInsert, particlesIndex_one_zero);
    EXPECT_TRUE(!d);

    vector<Particle*> particlesIndex_one_one = univ->getGrid()[1][1].getParticles();
    bool e  = checkParticleInCell(particleToInsert, particlesIndex_one_one);
    EXPECT_TRUE(!e);

    EXPECT_EQ(taillGrid(univ->getGrid()),4);
}





// test the changment of position : we put the particle in the in grid[1][1]
TEST_F(SetUpForAllTests, HandlesChangingPosition) {
    particleToInsert->setPosition(Vecteur(11,9,7));
    particleToInsert->setCellPositionChanged(true);
    updateGrid(*univ);


    vector<Particle*> particlesIndex_zero_zero = univ->getGrid()[0][0].getParticles();
    bool b  = checkParticleInCell(particleToInsert, particlesIndex_zero_zero);
    EXPECT_TRUE(!b);

    vector<Particle*> particlesIndex_zero_one = univ->getGrid()[0][1].getParticles();
    bool c  = checkParticleInCell(particleToInsert, particlesIndex_zero_one);
    EXPECT_TRUE(!c);

    vector<Particle*> particlesIndex_one_zero = univ->getGrid()[1][0].getParticles();
    bool d  = checkParticleInCell(particleToInsert, particlesIndex_one_zero);
    EXPECT_TRUE(d);

    vector<Particle*> particlesIndex_one_one = univ->getGrid()[1][1].getParticles();
    bool e  = checkParticleInCell(particleToInsert, particlesIndex_one_one);
    EXPECT_TRUE(!e);
}

// check the universe after the application of stromerVerletPotential function
TEST_F(SetUpForAllTests, HandlesStormerVerlet) {
    //ADD_FAILURE() << "universe debut : \n" << *univ ;
    stromerVerletPotential(*univ, 2,0.5, false, "");
    // The new expected position for the particle is grid[1][1]
    vector<Particle*> particlesIndex_zero_zero = univ->getGrid()[0][0].getParticles();
    bool b  = checkParticleInCell(particleToInsert, particlesIndex_zero_zero);
    EXPECT_TRUE(!b);

    vector<Particle*> particlesIndex_zero_one = univ->getGrid()[0][1].getParticles();
    bool c  = checkParticleInCell(particleToInsert, particlesIndex_zero_one);
    EXPECT_TRUE(!c);

    vector<Particle*> particlesIndex_one_zero = univ->getGrid()[1][0].getParticles();
    bool d  = checkParticleInCell(particleToInsert, particlesIndex_one_zero);
    EXPECT_TRUE(!d);

    vector<Particle*> particlesIndex_one_one = univ->getGrid()[1][1].getParticles();
    bool e  = checkParticleInCell(particleToInsert, particlesIndex_one_one);
    EXPECT_TRUE(e);
}

// check the universe after the application of stromerVerletPotential function
// when there is 2 particles
TEST_F(SetUpForAllTests, Handles2Particles) {
    secondParticle->setCellPositionChanged(true);
    updateGrid(*univ);
    ADD_FAILURE() << "universe debut : \n" << *univ ;
    /*
    stromerVerletPotential(*univ, 2,0.5, false, "");
    // The new expected position for the particle is grid[1][1]
    vector<Particle*> particlesIndex_zero_zero = univ->getGrid()[0][0].getParticles();
    bool b  = checkParticleInCell(particleToInsert, particlesIndex_zero_zero);
    EXPECT_TRUE(!b);

    vector<Particle*> particlesIndex_zero_one = univ->getGrid()[0][1].getParticles();
    bool c  = checkParticleInCell(particleToInsert, particlesIndex_zero_one);
    EXPECT_TRUE(!c);

    vector<Particle*> particlesIndex_one_zero = univ->getGrid()[1][0].getParticles();
    bool d  = checkParticleInCell(particleToInsert, particlesIndex_one_zero);
    EXPECT_TRUE(!d);

    vector<Particle*> particlesIndex_one_one = univ->getGrid()[1][1].getParticles();
    bool e  = checkParticleInCell(particleToInsert, particlesIndex_one_one);
    EXPECT_TRUE(e);
    */
}


