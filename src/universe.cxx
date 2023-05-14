#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "universe.hxx"

Universe::Universe(std::vector<Particle *> &particles, int dim, double rCut, double *LD) : particles(particles), rCut(rCut)
{
    // Organize the particles in cells
    this->dim = dim;
    // lD
    this->lD = new double[dim];
    for (int i = 0; i < dim; i++)
    {
        this->lD[i] = LD[i];
    }
    // nCD
    this->nCD = new long[dim];
    for (int i = 0; i < dim; i++)
    {
        this->nCD[i] = (long)(lD[i] / rCut);
    }
    // grid
    for (auto &particle : particles)
    {
        int x_cell = (int)(particle->getPosition().getX() / rCut);
        int y_cell = (int)(particle->getPosition().getY() / rCut);
        int z_cell = (int)(particle->getPosition().getZ() / rCut);
        // calculate indexs of a particle in the grid
        std::vector<int> indexs = {x_cell, y_cell, z_cell};
        // add the particle to the vector of particles corresponding to that indexs
        auto it = grid.find(indexs);
        if (it != grid.end())
        {
            it->second.push_back(particle);
        }
        else
        {
            std::vector<Particle *> particle_in_indexs = {particle};
            grid[indexs] = particle_in_indexs;
        }
    }
};

ostream &operator<<(std::ostream &os, const Universe &universe)
{
    for (const auto &grid_element : universe.grid)
    {
        os << "particles in cell ( ";
        const auto &grid_index = grid_element.first;
        auto it_indexs = grid_index.begin();
        while (it_indexs != std::prev(grid_index.end()))
        {
            os << *it_indexs << " , ";
            ++it_indexs;
        }
        os << *it_indexs << " ) are : ";
        os << "[ (";
        const auto &particles = grid_element.second;
        auto it_particles = particles.begin();
        while (it_particles != std::prev(particles.end()))
        {
            os << (*it_particles)->getPosition().getX() << " , " << (*it_particles)->getPosition().getY() << " , " << (*it_particles)->getPosition().getZ() << " ) , ( ";
            ++it_particles;
        }
        os << (*it_particles)->getPosition().getX() << " , " << (*it_particles)->getPosition().getY() << " , " << (*it_particles)->getPosition().getZ() << " ) ]" << std::endl;
    }
    return os;
}

/**
 * Calculates the potential interaction forces between a particle and its neighboring particles
 * @param universe
 */
void calculate_interaction_forces_potentiel(Universe universe)
{
    vector<Particle *> particles = universe.particles;
    double rCut = universe.rCut;
    std::map<std::vector<int>, std::vector<Particle *>> grid = universe.grid;

    std::vector<int> indexs_voisins = {-1, 0, 1};
    for (auto &particule : particles)
    {
        for (auto i : indexs_voisins)
        {
            for (auto j : indexs_voisins)
            {
                for (auto k : indexs_voisins)
                {
                    int x_cell = (int)(particule->getPosition().getX() / rCut);
                    int y_cell = (int)(particule->getPosition().getY() / rCut);
                    int z_cell = (int)(particule->getPosition().getZ() / rCut);
                    std::vector<int> indexs = {x_cell + i, y_cell + j, z_cell + k};
                    // add the particle to the vector of particles corresponding to that indexs
                    auto it = grid.find(indexs);
                    if (it != grid.end())
                    {
                        for (auto &particule_voisin : it->second)
                        {
                            if (particule != particule_voisin)
                            {
                                double rIJ = particule->getPosition().distanceToVect(particule_voisin->getPosition());
                                if (rIJ <= rCut)
                                {
                                    Vecteur fIJ1 = particule->getPosition();
                                    Vecteur fIJ2 = particule_voisin->getPosition();
                                    Vecteur fIJ = fIJ1 - fIJ2;
                                    fIJ = fIJ * 24 * universe.epsilon * (1 / pow(rIJ, 2)) * pow(universe.sigma / rIJ, 6) *
                                          (1 - (2 * pow(universe.sigma / rIJ, 6)));
                                    particule->force = particule->force + fIJ;
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
void modify_grid(Universe universe, Particle particleI, Vecteur ancien_position)
{
    double rCut = universe.rCut;
    for (auto &particle : universe.particles)
    {
        int x_cell = (int)(particle->getPosition().getX() / rCut);
        int y_cell = (int)(particle->getPosition().getY() / rCut);
        int z_cell = (int)(particle->getPosition().getZ() / rCut);
        std::vector<int> new_grid_position = {x_cell, y_cell, z_cell};
        std::vector<int> ancien_grid_position = {(int)(ancien_position.getX() / rCut), (int)(ancien_position.getY() / rCut), (int)(ancien_position.getZ() / rCut)};

        // we remove the particle from her previous cell
        auto it_ancien = universe.grid.find(ancien_grid_position);
        vector<Particle *> ancienParticleList = it_ancien->second;
        std::vector<Particle *>::iterator itera = ancienParticleList.begin();
        while (itera != ancienParticleList.end())
        {
            if (*itera == particle)
            {
                itera = ancienParticleList.erase(itera); // Deletes the element pointed to by it
            }
            else
            {
                ++itera;
            }
        }

        // we put the particle in the new cell
        auto it = universe.grid.find(new_grid_position);
        if (it != universe.grid.end())
        {
            it->second.push_back(particle);
        }
        else
        {
            std::vector<Particle *> particle_in_indexs = {particle};
            universe.grid[new_grid_position] = particle_in_indexs;
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

/**
 * Uses Stromer Verlet integration scheme to simulate the movement of particles with potential force
 * @param universe
 * @param tEnd
 * @param deltaT
 * @param outputStream
 */
void stromerVerlet_potentiel(Universe universe, double tEnd, double deltaT,
                             ofstream &outputStream)
{
    // Calculate the initial forces
    vector<Particle *> particleList = universe.particles;
    vector<Vecteur> fOld;
    calculate_interaction_forces_potentiel(universe);
    for (auto &particleI : particleList)
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

            Vecteur copy_ancien_position = particleI->position;
            particleI->position =
                particleI->position + (particleI->vitesse +
                                       particleI->force * (0.5 / particleI->masse) * deltaT) *
                                          deltaT;
            modify_grid(universe, *particleI, copy_ancien_position);
            fOld[index] = particleI->force;
            index++;
        }

        calculate_interaction_forces_potentiel(universe);

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
