#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "universe.hxx"
#include "cell.hxx"

Universe::Universe(std::vector<Particle *> &particles, int dim, double rCut, double *lD) : particles(particles), rCut(rCut)
{
    // Organize the particles in cells
    this->dim = dim;
    // lD
    this->lD = new double[dim];
    for (int i = 0; i < dim; i++)
    {
        this->lD[i] = lD[i];
    }
    // nCD
    this->nCD = new long[dim];
    for (int i = 0; i < dim; i++)
    {
        this->nCD[i] = (long)(lD[i] / rCut);
    }
    //grid
    //Init the grid to the proper size x = nCD[0], y = nCD[1], z = nCD[2]
    this->grid.resize(nCD[0]);
    for (auto &yAxe : grid){
        yAxe.resize(nCD[1]);
    }
    //Filling the grid
    for (auto particle : particles)
    {
        int x_cell = (int)(particle->getPosition().getX() / rCut);
        int y_cell = (int)(particle->getPosition().getY() / rCut);
        int z_cell = (int)(particle->getPosition().getZ() / rCut);
        grid[x_cell][y_cell].particles.push_back(particle);//2D
    }
};



ostream &operator<<(std::ostream &os, const Universe &universe)
{
    for (int x = 0; x < universe.nCD[0]; x ++){
        for (int y = 0; y < universe.nCD[1]; y ++){
            os << "Cell:" << x << "," << y << "particles";
            for (auto particle : universe.grid[x][y].getParticles()){
                os << particle << " ; ";
            }
            os << endl;
        }
    }

    return os;
}


/**
 * Calculates the potential interaction forces between a particle and its neighboring particles
 * @param universe
 */
void interactionForcesPotentiel(Universe & universe)
{
    for (int x = 0; x < universe.nCD[0]; x++){
        for (int y = 0; y < universe.nCD[1]; y++){
            //We are in the cell (x,y)
            //We loop over the neighboring cells, if they exist we consider the particles inside
            for (auto particleI : universe.grid[x][y].getParticles()){
                for (int dx = -1; dx <= 1; dx ++){
                    for (int dy = -1; dy <= 1; dy++){
                        if (x + dx < universe.nCD[0] && x + dx >= 0 && y + dy < universe.nCD[1] && y + dy >= 0){
                            for (auto particleJ : universe.grid[x][y].getParticles()){
                                if (!(particleI == particleJ)){
                                    Vecteur rIJVect = particleJ->position - particleI->position;
                                    double rIJ = rIJVect.norm();
                                    if (rIJ < universe.rCut && rIJ > 0){
                                        double fIJMagnitude = 24 * universe.epsilon * (1/(rIJ*rIJ)) *
                                            pow((universe.sigma /rIJ), 6) * (1- 2 * pow((universe.sigma/rIJ), 6));
                                        particleI->force = particleI->force + rIJVect * fIJMagnitude;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


/**
 * Called after change of position of particles, it removes a particle from the old cell and puts it in the new cell
 * @param universe
 * @param particleI
 * @param ancien_position
 */
void updateGrid(Universe &universe)
{
    for (int x = 0; x < universe.nCD[0]; x++){
        for (int y = 0; y < universe.nCD[1]; y++){
            for (auto particle : universe.grid[x][y].particles){
                if (particle->isCellPositionChanged()){
                    int x_cell = (int)(particle->getPosition().getX() / universe.rCut);
                    int y_cell = (int)(particle->getPosition().getY() / universe.rCut);
                    int z_cell = (int)(particle->getPosition().getZ() / universe.rCut);
                    //2D
                    universe.grid[x_cell][y_cell].removeParticle(particle);
                    universe.grid[x_cell][y_cell].particles.push_back(particle);//2D
                    //We set back the position to its initial state
                    particle->setCellPositionChanged(false);
                }
            }
        }
    }
}

/**
 * Uses Stromer Verlet integration scheme to simulate the movement of particles with potential force
 * @param universe
 * @param tEnd
 * @param deltaT
 * @param outputStream
 */
void stromerVerletPotential(Universe &universe, double tEnd, double deltaT,
                            ofstream &outputStream)
{
    // Calculate the initial forces
    vector<Particle *> particleList = universe.particles;
    vector<Vecteur> fOld;
    interactionForcesPotentiel(universe);
    for (auto particleI : particleList)
    {
        fOld.push_back(particleI->force);
    }

    double t = 0.0;
    while (t < tEnd)
    {
        t += deltaT;
        int index = 0;

        for (auto &particleI : particleList)
        {
            particleI->position =
                particleI->position + (particleI->vitesse +
                                       particleI->force * (0.5 / particleI->masse) * deltaT) *
                                      deltaT;
            //The position is changed -> declare the change
            particleI->setCellPositionChanged(true);
            updateGrid(universe);

            fOld[index] = particleI->force;
            index++;
        }

        interactionForcesPotentiel(universe);

        index = 0; // will be reused after this
        for (auto &particleI : particleList)
        {
            particleI->vitesse = particleI->vitesse + (particleI->force + fOld[index]) *
                                                      deltaT * (0.5 / particleI->masse);
            index++;
        }
        for (auto &particleI : particleList)
        {
            outputStream << particleI->position << endl;
        }
    }
}

void Universe::calculateForcesUni()
{
    calculateForces(particles);
}

void Universe::calculateForcesSlowUni()
{
    calculateForcesSlow(particles);
}




