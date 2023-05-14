#ifndef CELL_H
#define CELL_H

#include <vector>

#include "particule.hxx"

class Cell {
private :
        vector<Particle *> particles;
        vector<Cell*> neighborCells;
        vector<long> position; //Position in the cells grid, pointer to array
public :
        Cell(vector<Particle *> &particles, vector<long> position): particles(particles), position(position){
                neighborCells.push_back(this);//The neighborCells contains the cell itself
        }
        vector<long> getPosition(){return position;}
        void setPosition(vector<long> position){this->position = position;}

        void addNeighborCell(Cell * cell){
                neighborCells.push_back(cell);
        };
};

#endif