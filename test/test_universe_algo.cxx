#include <iostream>
#include <random>
#include <list>

#include "universe.hxx"
#include "particule.hxx"


using namespace std;


vector<Particle *> constructUniverse(double k){
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);

    vector<Particle*> particleList;

    //Make sure what's the pow function is encodable
    auto start_time = chrono::steady_clock::now();
    for (int i = 0; i < pow(pow(2.0,k), 3.0); ++i){
        string type = "(nil)";
        //Using auto deduces automatically the type
        Particle * particleToInsert = new Particle(Vecteur(dist(mt), dist(mt), dist(mt)), Vecteur(7.0, 7.0, 7.0)
                , 1.0, i, Vecteur(0.0, 0.0, 0.0), type);
        particleList.push_back(particleToInsert);
    }
    auto exec_time = chrono::steady_clock::now() - start_time;

    cout << "Time taken for the insertion with size = 2**(3 * " << k << ") ----> "
         << chrono::duration<double, milli>(exec_time).count() << "ms" << endl;

    return particleList;
}


int main() {

    //TP3-Question7: test insertion en utilisant les listes


    // we will test class Universe and see the influence of rCut on the distribution of particles in the grid
    double rCut = 0.1; //Cutoff distance for potential
    //int dim = 3 ; // dimension
    //double lD[3] = {10,10,10}; //Characteristic length depending on direction
    int dim = 1; // dimension
    double lD[1] = {10}; //Characteristic length depending on direction
    vector<Particle *> particles2 = constructUniverse(2);
    Universe *univ2 = new Universe(particles2, dim, 0.3, lD);

    ofstream outputFile;
    outputFile.open("../test_universe_algo.txt");
    stromerVerlet_potentiel(*univ2, 468.5, 0.015, outputFile);
    outputFile.close();

}