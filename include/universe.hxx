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


class Universe {
private :
        double epsilon;
        double sigma;
        std::vector<Particle *> particles;
        double rCut; //Cutoff distance for potential
        int dim;
        double *lD; //Characteristic length depending on direction
        long *nCD;
//        std::map<std::vector<int>, std::vector<Particle *>> grid;
        std::vector<std::vector<Cell>> grid;
        //TODO: Use templates to define these
public:
        Universe(std::vector<Particle*>& particles, int dim, double rCut, double* LD);

//        friend void calculate_interaction_forces_potentiel(Universe universe);
//        friend void stromerVerlet_potentiel(Universe universe, double tEnd, double deltaT,
//                                 ofstream outputStream);
//        friend void modify_grid(Universe universe, Particle particleI, Vecteur ancien_position) ;
//
//        friend void calculateForces(vector<Particle *> &);
        friend ostream& operator<<(ostream &o, const Universe &);
//
//        const vector<Particle *> &getParticles() const;

        friend void interactionForcesPotentiel(Universe & universe);
        friend void stromerVerletPotential(Universe & universe, double tEnd, double deltaT, ofstream &outputStream);
        friend void updateGrid(Universe & universe);

        void calculateForcesUni();
        void calculateForcesSlowUni();
};

const vector<Particle *> &Cell::getParticles() const {
    return particles;
}

ostream& operator<<(ostream & o, const Universe & u);





#endif