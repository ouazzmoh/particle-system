#include <iostream>
#include <random>
#include <list>

#include "universe.hxx"
#include "particule.hxx"




int main() {

    double lD[2] = {250, 100};
    double epsilon = 5;
    double sigma = 1;
    double rCut = 2.5*sigma;
    double d_t = 0.00005;

    vector<Particle * > particleList;
    double spacing = pow(2, 1/6) /sigma;

    int id = 0;
    //Red square
    double xInit = 95;
    double yInit = 45 + spacing;

    for (int i = 0; i < 40; i++){
        for ( int j = 0; j < 40; j ++) {
            string type = "nil";
                Particle *particleToInsert = new Particle(Vecteur(xInit + i * spacing, yInit + j * spacing, 0.0),
                                                          Vecteur(0.0, -10.0, 0.0), 1.0, id, Vecteur(0.0, 0.0, 0.0),
                                                          type);
                particleList.push_back(particleToInsert);
            id ++;
        }
    }
    //Blue Rectangle
    double xInit2 = 35;
    double yInit2 = 5;
    for (int i = 0; i < 160; i++){
        for (int j = 0; j < 40; j ++) {
            string type = "nil";
            Particle *particleToInsert = new Particle(Vecteur(xInit2 + i * spacing, yInit2 + j * spacing, 0.0),
                                                      Vecteur(0.0, 0.0, 0.0), 1.0, id, Vecteur(0.0, 0.0, 0.0),
                                                      type);
            particleList.push_back(particleToInsert);
            id++;
        }

    }



//    ofstream init;
//    init.open("../../demo/tp4_init.vtu");
//    printVtk(particleList, init);
//    init.close();

    Universe *univ = new Universe(particleList, 2, rCut, lD, epsilon, sigma);



    stromerVerletPotential(*univ, 19.5, d_t, true, "../../demo/tp4_application/", false, 0, 0);


    return 0;
};