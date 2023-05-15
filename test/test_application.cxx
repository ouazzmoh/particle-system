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
    double d_t = 0.5;

    vector<Particle * > particleList;

    //Blue Square
    for (int i = 0; i < 160; i++){
        for ( int j = 0; j < 40; j ++) {
            string type = "nil";
            Particle *particleToInsert = new Particle(Vecteur(i , j, 0.0),
                                                      Vecteur(0.0, 0.0, 0.0), 1.0, i, Vecteur(0.0, 0.0, 0.0), type);
            particleList.push_back(particleToInsert);
        }

    }
    //Red square
    for (int i = 0; i < 40; i++){
        for ( int j = 0; j < 40; j ++) {
            string type = "nil";
            Particle *particleToInsert = new Particle(Vecteur(60 + i , pow(2, 1/6)/sigma + 40 + j, 0.0),
                                                      Vecteur(0.0, 0.0, 0.0), 1.0, i, Vecteur(0.0, 0.0, 0.0), type);
            particleList.push_back(particleToInsert);
        }
    }


    Universe *univ = new Universe(particleList, 2, rCut, lD);

    ofstream outputFile;
    outputFile.open("../../test_application.txt");
     for (auto &particle: particleList){
         outputFile << particle->getPosition() << endl;
     }
    stromerVerletPotential(*univ, 19.5, d_t, outputFile);
    outputFile.close();

    return 0;
};