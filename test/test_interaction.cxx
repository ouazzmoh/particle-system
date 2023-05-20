#include "universe.hxx"
#include "particule.hxx"




int main() {
    double lD[2] = {100, 100};
    double epsilon = 5;
    double sigma = 1;
    double rCut = 10;

    vector<Particle * > particleListAbs;

    string type = "nil";
    //Absorption
    Particle *particleToInsert1 = new Particle(Vecteur(50.0 - 10, 50.0, 0.0),
                                              Vecteur(5.0, 0.0, 0.0), 1.0, 1, Vecteur(0.0, 0.0, 0.0), type);
    Particle *particleToInsert2 = new Particle(Vecteur(50.0 + 10, 50.0, 0.0),
                                               Vecteur(-5.0, 0.0, 0.0), 1.0, 2, Vecteur(0.0, 0.0, 0.0), type);
    particleListAbs.push_back(particleToInsert1);
    particleListAbs.push_back(particleToInsert2);
    Universe *univAbs = new Universe(particleListAbs, 2, rCut, lD, epsilon, sigma, 0);
    startSimulation(*univAbs, 20, 0.00005, true, "../../demo/interSimple/", true, 0, 0);

    return 0;
};



