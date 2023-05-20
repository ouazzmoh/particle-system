#ifndef CELL_H
#define CELL_H

#include <unordered_set>

#include "particule.hxx"

class Cell {
private :
        unordered_set<Particle * > particles;
public :

        //Default constructor, useful to create the grid and initializes to empty vector
        Cell() : particles({}) {}

        Cell(unordered_set<Particle *> &particles): particles(particles){}

        void removeParticle(Particle * particle);

        friend class Universe;
        friend void updateGrid(Universe &universe, Particle *);
        friend void updateGrid1D(Universe &universe, Particle * particle);
        friend void updateGrid3D(Universe &universe, Particle * particle);

        const unordered_set<Particle *> &getParticles() const;

};

#endif