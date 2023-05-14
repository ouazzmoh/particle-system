//
// Created by Simo Ouazzani  on 13/04/2023.
//

#include <iostream>
#include <list>
#include <random>

#include "universe.hxx"


int main(){

        //TP4:Question1 systeme a deux particules on trace le potentiel en fonction de la distance

        ofstream outputQ1;
        outputQ1.open("../plot2particleSystem.txt");//By default it puts it in build/test folder
        double sigma = 1.0;
        double epsilon = 1.0;
        double r = 0.001;
        while (r < 0.3){
                outputQ1 << r << "      " << 4 * epsilon* pow((sigma/r), 6) * (pow((sigma/r), 6) - 1) << endl;
                r += 0.001;
        }
        outputQ1.close();
        //Plot with gnuplot : voir uFCTr.png

        //Nous pourrons considerer que le potentiel est nul a partir de r = 0.15


        vector<Particle * > particleList;
        random_device rd;
        mt19937 mt(rd());
        uniform_real_distribution<double> dist(0.0, 100.0);
        for (int i = 0; i < 1000; i++){
                string type = "(nil)";
                Particle* particleToInsert = new Particle(Vecteur(dist(mt), dist(mt), dist(mt)), Vecteur(0.0, 0.0, 0.0)
                        , 1.0, i, Vecteur(0.0, 0.0, 0.0), type);

                particleList.push_back(particleToInsert);
        }



        //Freeing all dynamically allocated memory
        for (auto particle: particleList){
                delete particle;
        }
        return 0;
}
