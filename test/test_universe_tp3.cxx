#include <iostream>
#include <random>
#include <list>

#include "universe.hxx"

using namespace std;

/**
 * TP3- Question 7
 * Construct a universe of size (2**k)**3, prints the insertion execution results
 * @param k : a double, for and int 2 cast to double before so it is encodable
 * @return Pointer to the universe, dynamic allocation
 *
 */
Universe * constructUniverse(double k){
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);
    double lD[2] = {10, 10};

    vector<Particle*> particleList;

    //Make sure what's the pow function is encodable
    auto start_time = chrono::steady_clock::now();
    for (int i = 0; i < fastPow(fastPow(2.0,k), 3.0); ++i){
        string type = "(nil)";
        //Using auto deduces automatically the type
        Particle * particleToInsert = new Particle(Vecteur(dist(mt), dist(mt), dist(mt)), Vecteur(7.0, 7.0, 7.0)
            , 1.0, i, Vecteur(0.0, 0.0, 0.0), type);
        particleList.push_back(particleToInsert);
    }
    auto exec_time = chrono::steady_clock::now() - start_time;

    cout << "Time taken for the insertion with size = 2**(3 * " << k << ") ----> "
         << chrono::duration<double, milli>(exec_time).count() << "ms" << endl;

    return new Universe(particleList, 2, 1, lD, 0, 0);
}


int main() {


    //TP3-Question 8 - 9 : Calcule des interactions
    //Pour diviser le temps de calcul en 2, On utilise la 3eme loi de Newton
    //en calculant les forces d'interactions symetriques, cad Fij = - Fji,
    for (int k = 3 ; k < 6; k++){
        cout << "-------------- k = " << k << "-------------------" << endl;
        Universe * universeK = constructUniverse(k);
        //TP3-Question 8 - 9 : Calcule des interactions
        //Pour diviser le temps de calcul en 2, On utilise la 3eme loi de Newton
        //en calculant les forces d'interactions, cad Fij = - Fji,
        auto start_time = chrono::steady_clock::now();
        universeK->calculateForcesSlowUni();
        auto exec_time = chrono::steady_clock::now() - start_time;
        cout << "Time taken for forces SLOW calculation with size = 2**(3 * " << k << ") ----> "
             << chrono::duration<double, milli>(exec_time).count() << "ms" << endl;
        auto start_time2 = chrono::steady_clock::now();
        universeK->calculateForcesUni();
        auto exec_time2 = chrono::steady_clock::now() - start_time2;
        cout << "Time taken for forces FASTER calculation with size = 2**(3 * " << k << ") ----> "
             << chrono::duration<double, milli>(exec_time2).count() << "ms" << endl;
    }

    //Question 10 :
    //TODO: Think about more ideas to optimize
    //Parmi les methodes pour optimiser c'est de calculer les forces que pour les particules qui sont proches.
    //Un algorithme possible: (composantes connexes)
    //      *Definir une distance minimal pour laquelle on associe les particules
    //      *Ne pas calculer les forces que a l'interieur de ces composantes connexes
    //      TODO: Implement this idea ??

}