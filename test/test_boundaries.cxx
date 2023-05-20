#include "universe.hxx"
#include "particule.hxx"


using namespace std;

int main() {
    /**
     * This test creates universes, and tests different boundary conditions on them
     */

    double lD[2] = {100, 100};
    double epsilon = 5;
    double sigma = 1;
    double rCut = 10;

    vector<Particle * > particleListAbs;
    vector<Particle * > particleListRef; //Reflexion
    vector<Particle * > particleListPer; //Periodic
    vector<Particle * > particleListRefLJ; //Reflexion using the LJ potential

    string type = "nil";
    //Absorption
    for (int i = 0; i < 4; i++){
        Particle *particleToInsert = new Particle(Vecteur(50.0 + 10 * i, 50.0, 0.0),
                                                  Vecteur(5.0, -4.0, 0.0), 1.0, 1, Vecteur(0.0, 0.0, 0.0), type);
        particleListAbs.push_back(particleToInsert);
    }
    Universe *univAbs = new Universe(particleListAbs, 2, rCut, lD, epsilon, sigma, 0);
    stromerVerletPotential(*univAbs, 20, 0.005, true, "../../demo/abs/", false, 0, 0);

    //Reflexion
    for (int i = 0; i < 4; i++){
        Particle *particleToInsert = new Particle(Vecteur(50.0 + 10 * i, 50.0, 0.0),
                                                  Vecteur(5.0, -4.0, 0.0), 1.0, 1, Vecteur(0.0, 0.0, 0.0), type);
        particleListRef.push_back(particleToInsert);
    }
    Universe *univRef = new Universe(particleListRef, 2, rCut, lD, epsilon, sigma, -1);
    stromerVerletPotential(*univRef, 20, 0.005, true, "../../demo/ref/", false, 0, 0);

    //Periodic
    for (int i = 0; i < 4; i++){
        Particle *particleToInsert = new Particle(Vecteur(50.0 + 10 * i, 50.0, 0.0),
                                                  Vecteur(5.0, -4.0, 0.0), 1.0, 1, Vecteur(0.0, 0.0, 0.0), type);
        particleListPer.push_back(particleToInsert);
    }
    Universe *univPer = new Universe(particleListPer, 2, rCut, lD, epsilon, sigma, 1);
    stromerVerletPotential(*univPer, 20, 0.005, true, "../../demo/per/", false, 0, 0);

    //RefLJlexion
    for (int i = 0; i < 4; i++){
        Particle *particleToInsert = new Particle(Vecteur(50.0 + 10 * i, 50.0, 0.0),
                                                  Vecteur(5.0, -5.0, 0.0), 1.0, 1, Vecteur(0.0, 0.0, 0.0), type);
        particleListRefLJ.push_back(particleToInsert);
    }
    Universe *univRefLJ = new Universe(particleListRefLJ, 2, rCut, lD, epsilon, sigma, 0);
    stromerVerletPotential(*univRefLJ, 20, 0.00005, true, "../../demo/refLj/", true, 0, 0);


    delete univAbs;
    delete univPer;
    delete univRef;
    delete univRefLJ;
    for (auto particle : particleListAbs){
        delete particle;
    }
    for (auto particle : particleListRef){
        delete particle;
    }
    for (auto particle : particleListRefLJ){
        delete particle;
    }
    for (auto particle : particleListPer){
        delete particle;
    }



    return 0;
};



