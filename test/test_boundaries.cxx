#include "universe.hxx"
#include "particule.hxx"




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
    stromerVerletPotential(*univAbs, 20, 0.005, true, "../../demo/abs/", false, false);

    //Reflexion
    for (int i = 0; i < 4; i++){
        Particle *particleToInsert = new Particle(Vecteur(50.0 + 10 * i, 50.0, 0.0),
                                                  Vecteur(5.0, -4.0, 0.0), 1.0, 1, Vecteur(0.0, 0.0, 0.0), type);
        particleListRef.push_back(particleToInsert);
    }
    Universe *univRef = new Universe(particleListRef, 2, rCut, lD, epsilon, sigma, -1);
    stromerVerletPotential(*univRef, 20, 0.005, true, "../../demo/ref/", false, false);

    //Periodic
    for (int i = 0; i < 4; i++){
        Particle *particleToInsert = new Particle(Vecteur(50.0 + 10 * i, 50.0, 0.0),
                                                  Vecteur(5.0, -4.0, 0.0), 1.0, 1, Vecteur(0.0, 0.0, 0.0), type);
        particleListPer.push_back(particleToInsert);
    }
    Universe *univPer = new Universe(particleListPer, 2, rCut, lD, epsilon, sigma, 1);
    stromerVerletPotential(*univPer, 20, 0.005, true, "../../demo/per/", false, false);

    //RefLJlexion
    for (int i = 0; i < 4; i++){
        Particle *particleToInsert = new Particle(Vecteur(50.0 + 10 * i, 50.0, 0.0),
                                                  Vecteur(5.0, -4.0, 0.0), 1.0, 1, Vecteur(0.0, 0.0, 0.0), type);
        particleListRefLJ.push_back(particleToInsert);
    }
    Universe *univRefLJ = new Universe(particleListRefLJ, 2, rCut, lD, epsilon, sigma, 0);
    stromerVerletPotential(*univRefLJ, 20, 0.00005, true, "../../demo/refLj/", true, false);

    return 0;
};



