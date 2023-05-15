#ifndef CELL_H
#define CELL_H

#include <vector>

#include "particule.hxx"

class Cell {
private :
        vector<Particle *> particles;
public :

        //Default constructor, useful to create the grid and initializes to empty vector
        Cell() : particles({}) {}

        Cell(vector<Particle *> &particles): particles(particles){}

        friend class Universe;
        friend void updateGrid(Universe &universe);

        const vector<Particle *> &getParticles() const;

        void removeParticle(Particle * particle){
            particles.erase(std::remove(particles.begin(), particles.end(), particle), particles.end());
        }
};

#endif