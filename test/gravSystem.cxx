#include <iostream>


#include "vecteur.hxx"
#include "particule.hxx"


using namespace std;


int main() {

        //Moving the particles
        vector<Particle *> gravitationalSystem;

        auto soleil = new Particle(Vecteur(0.0, 0.0, 0.0),
                                  Vecteur(0.0, 0.0, 0.0), 1.0, 0,
                                  Vecteur(0.0, 0.0, 0.0), "Soleil");
        auto terre = new Particle(Vecteur(0.0, 1.0, 0.0),
                                  Vecteur(-1.0, 0.0, 0.0),3e-6, 1,
                                  Vecteur(0.0, 0.0, 0.0), "Terre");
        auto jupiter = new Particle(Vecteur(0.0, 5.36, 0.0),
                                    Vecteur(-0.425, 0.0, 0.0), 9.55e-4, 2,
                                    Vecteur(0.0, 0.0, 0.0), "Jupiter");
        auto haley = new Particle(Vecteur(34.75, 0.0, 0.0),
                                  Vecteur(0.0, 0.0296, 0.0), 1e-14, 3,
                                  Vecteur(0.0, 0.0, 0.0), "Haley");



        gravitationalSystem.insert(gravitationalSystem.end(), soleil);
        gravitationalSystem.insert(gravitationalSystem.end(), terre);
        gravitationalSystem.insert(gravitationalSystem.end(), jupiter);
        gravitationalSystem.insert(gravitationalSystem.end(), haley);


        ofstream outputFile;
        outputFile.open("../../coords.txt");
        stromerVerlet(gravitationalSystem, 468.5, 0.015, outputFile);
        outputFile.close();

        //To plot the result we use gnuplot
        //Check gravSystem.png for results

        delete soleil;
        delete terre;
        delete jupiter;
        delete haley;
}