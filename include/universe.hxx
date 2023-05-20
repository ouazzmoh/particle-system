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
        std::vector<Cell> grid1D;
        std::vector<std::vector<Cell>> grid;
        std::vector<std::vector<std::vector<Cell>>> grid3D;
public:
        void setGrid1D(const std::vector<Cell> &grid1D);

        void setGrid3D(const std::vector<std::vector<std::vector<Cell>>> &grid3D);

private:
        int boundCond = 0; //boundary condition, 0:absorption/ 1:periodic (portal) / -1:reflection
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
        friend std::ostream& operator<<(std::ostream &o, const Universe &);


        friend void interactionForcesPotentiel(Universe & universe, bool ljReflexion, double G);
        friend void interactionForcesPotentiel3D(Universe & universe, bool ljReflexion, double G);
        friend void interactionForcesPotentiel1D(Universe & universe, bool ljReflexion, double G);
        friend void stromerVerletPotential(Universe & universe, double tEnd, double deltaT,bool visual, std::string path,
                                           bool ljReflexion, double G, double eCD);
        friend void stromerVerletPotential3D(Universe & universe, double tEnd, double deltaT,bool visual,  std::string path,
                                             bool ljReflexion, double G, double eCD);
        friend void stromerVerletPotential1D(Universe & universe, double tEnd, double deltaT,bool visual,  std::string path,
                                             bool ljReflexion, double G, double eCD);
        friend void startSimulation(Universe & universe, double tEnd, double deltaT,bool visual,  std::string path,
                                             bool ljReflexion, double G, double eCD);
        friend void updateGrid(Universe & universe, Particle *);
        friend void updateGrid1D(Universe &universe, Particle * particle);
        friend void updateGrid3D(Universe &universe, Particle * particle);
        void calculateForcesUni();
        void calculateForcesSlowUni();

        const std::vector<Particle *> &getParticles() const;
};



std::ostream& operator<<(std::ostream & o, const Universe & u);

void printVtk(std::vector<Particle *> particleList, std::ostream & outputStream);
void stromerVerletPotential(Universe & universe, double tEnd, double deltaT, bool visual,  std::string path,
                            bool ljReflexion, double G, double eCD);
void stromerVerletPotential3D(Universe & universe, double tEnd, double deltaT,bool visual,  std::string path,
                                     bool ljReflexion, double G, double eCD);
void stromerVerletPotential1D(Universe & universe, double tEnd, double deltaT,bool visual,  std::string path,
                                     bool ljReflexion, double G, double eCD);



void startSimulation(Universe & universe, double tEnd, double deltaT,bool visual,  std::string path,
                              bool ljReflexion, double G, double eCD);



#endif