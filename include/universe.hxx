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
        std::vector<std::vector<Cell>> grid;
        int boundCond = 0; //boundary condition, 0:absorption/ 1:periodic (portal) / -1:reflection
        //TODO: Use templates to define these
public:
        Universe(std::vector<Particle*>& particles, int dim, double rCut, double* LD, double epsilon, double sigma, int boundCond =0);

        friend ostream& operator<<(ostream &o, const Universe &);

        friend void interactionForcesPotentiel(Universe & universe);
        friend void stromerVerletPotential(Universe & universe, double tEnd, double deltaT, ofstream &outputStream, string path);
        friend void updateGrid(Universe & universe);

        void calculateForcesUni();
        void calculateForcesSlowUni();

        const vector<Particle *> &getParticles() const;
};



ostream& operator<<(ostream & o, const Universe & u);

void printVtk(vector<Particle *> particleList, ostream & outputStream);
void stromerVerletPotential(Universe & universe, double tEnd, double deltaT, ofstream &outputStream, string path);

void updateGrid(Universe & universe);

#endif