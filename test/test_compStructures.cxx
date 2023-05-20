#include <iostream>

#include "vecteur.hxx"
#include "particule.hxx"

using namespace std;


int main() {
        list<int> particleTest = {64, 128, 1024, 2048};
        cout << "------------Time results for lists ------------" << endl;
        for (auto &testSize: particleTest) {
                auto start_time = chrono::steady_clock::now();
                constructParticleList(testSize, false);
//printing execution time
                auto exec_time = chrono::steady_clock::now() - start_time;
                cout << "Time taken for the test size of " << testSize << " :"
                     << chrono::duration<double, milli>(exec_time).count() << "ms" << endl;
        }
        cout << "------------Time results for sets ------------" << endl;
        for (auto &testSize: particleTest) {
                auto start_time = chrono::steady_clock::now();
                constructParticleSet(testSize, false);
//printing execution time
                auto exec_time = chrono::steady_clock::now() - start_time;
                cout << "Time taken for the test size of " << testSize << " :"
                     << chrono::duration<double, milli>(exec_time).count() << "ms" << endl;
        }


        cout << "------------Time results for vectors ------------" << endl;
        for (auto &testSize: particleTest) {
                auto start_time = chrono::steady_clock::now();
                constructParticleVector(testSize, false);
//printing execution time
                auto exec_time = chrono::steady_clock::now() - start_time;
                cout << "Time taken for the test size of " << testSize << " :"
                     << chrono::duration<double, milli>(exec_time).count() << "ms" << endl;
        }
        
        cout << "------------Time results for deques ------------" << endl;
        for (auto &testSize: particleTest) {
                auto start_time = chrono::steady_clock::now();
                constructParticleDeque(testSize, false);
//printing execution time
                auto exec_time = chrono::steady_clock::now() - start_time;
                cout << "Time taken for the test size of " << testSize << " :"
                     << chrono::duration<double, milli>(exec_time).count() << "ms" << endl;
        }

/**
 * Les résultats des listes sont meilleurs, car ce sont des listes doublement chainées avec un cout d'insertion o(1)
 * contrairement au set qui son des BST, donc cout d'inserion est de cout o(log(n))
 * Dans la suite nous travaillerons avec la structure list
 *
 * Nous utiliserons les vecteurs car pour des particules de grandes tailles, il nous donne de bonnes resultats
 */
}