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

        /**
         * Constructor for the universe
         * @param particles
         * @param dim
         * @param rCut
         * @param LD
         * @param epsilon
         * @param sigma
         * @param boundCond : by default the property of the universe is absorption
         */
        Universe(std::vector<Particle*>& particles, int dim, double rCut, double* LD, double epsilon, double sigma, int boundCond =0);
        std::vector<std::vector<Cell>> getGrid();
        friend ostream& operator<<(ostream &o, const Universe &);


        friend void interactionForcesPotentiel(Universe & universe, bool ljReflexion, double G);
        friend void stromerVerletPotential(Universe & universe, double tEnd, double deltaT,bool visual, string path,
                                           bool ljReflexion, double G, double eCD);
        friend void updateGrid(Universe & universe, Particle *);
        void calculateForcesUni();
        void calculateForcesSlowUni();

        const vector<Particle *> &getParticles() const;
};



ostream& operator<<(ostream & o, const Universe & u);

void printVtk(vector<Particle *> particleList, ostream & outputStream);
void stromerVerletPotential(Universe & universe, double tEnd, double deltaT, bool visual,  string path,
                            bool ljReflexion, double G, double eCD);



#endif