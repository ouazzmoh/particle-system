#ifndef UNIVERSE_H
#define UNIVERSE_H

#include <string>
#include <fstream>
#include <iostream>
#include <vector>


#include "particule.hxx"
#include "cell.hxx"
#include "vecteur.hxx"
#include <map>


using namespace std;

class Universe {
private :
        double epsilon;
        double sigma;
        std::vector<Particle *> particles;
        double rCut; //Cutoff distance for potential
        int dim;
        double *lD; //Characteristic length depending on direction
        long *nCD;
        std::map<std::vector<int>, std::vector<Particle *>> grid;
        //TODO: Use templates to define these
public:
        Universe(std::vector<Particle*>& particles, int dim, double rCut, double* LD);

        friend void calculate_interaction_forces_potentiel(Universe universe);
        friend void stromerVerlet_potentiel(Universe universe, double tEnd, double deltaT,
                                 ofstream &outputStream);
        friend void modify_grid(Universe universe, Particle particleI, Vecteur ancien_position) ;

        friend void calculateForces(vector<Particle *> &);
        friend ostream& operator<<(ostream &o, const Universe &);

        const vector<Particle *> &getParticles() const;

        void calculateForcesUni();
        void calculateForcesSlowUni();
};

ostream& operator<<(ostream & o, const Universe & u);





#endif