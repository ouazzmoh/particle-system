#include <iostream>
#include <random>
#include <list>

#include "universe.hxx"
#include "particule.hxx"




int main() {

    double lD[2] = {250, 40};
    double epsilon = 5;
    double sigma = 1;
    double rCut = 2.5*sigma;
    double d_t = 0.00005;

    vector<Particle * > particleList;

    //Red square
//    for (int i = 0; i < 10; i++){
//        for ( int j = 0; j < 10; j ++) {
    string type = "nil";
    Particle *particleToInsert = new Particle(Vecteur(60 , pow(2, 1/6)/sigma + 40, 0.0),
                                              Vecteur(0.0, 10.0, 0.0), 1.0, 1, Vecteur(0.0, 0.0, 0.0), type);
    particleList.push_back(particleToInsert);
//        }
//    }



    Universe *univ = new Universe(particleList, 2, rCut, lD, epsilon, sigma);

    ofstream init;
    ofstream end;
    init.open("../../demo/sim0.vtu");
    end.open("../../demo/sim1.vtu");

    printVtk(univ->getParticles(), init);

    stromerVerletPotential(*univ, 1, d_t, init, "../../demo/");

    printVtk(univ->getParticles(), end);
    init.close();
    end.close();

    return 0;
};