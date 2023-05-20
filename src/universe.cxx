#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

#include "universe.hxx"
#include "cell.hxx"

using namespace std;


Universe::Universe(std::vector<Particle *> &particles, int dim, double rCut, double *lD, double epsilon, double sigma, int boundCond) : particles(particles), rCut(rCut), epsilon(epsilon),
sigma(sigma), boundCond(boundCond)
{
    if (rCut <= 0){
        throw invalid_argument("the value of rCut needs to be positive");
    }

    if (!(boundCond == 1 || boundCond == 0 ||boundCond == -1)){
        throw invalid_argument("the value of boundCond needs to be either 0(Absorption), 1(Periodic), -1(Reflexive)");
    }

    // Organize the particles in cells
    this->dim = dim;
    // lD
    this->lD = new double[dim];
    for (int i = 0; i < dim; i++)
    {
        if (lD[i] > 0 && lD[i] >= rCut)
        {
            this->lD[i] = lD[i];
        }
        else{
            throw invalid_argument("the value of lD needs to be positive and superior to rCut");
        }
    }
    // nCD
    this->nCD = new long[dim];
    for (int i = 0; i < dim; i++)
    {
        this->nCD[i] = (long)(lD[i] / rCut);
    }
    //grid
    if (dim == 2) {
        //grid
        //Init the grid to the proper size x = nCD[0], y = nCD[1], z = nCD[2]
        this->grid.resize(nCD[0]);
        for (auto &yAxe: grid) {
            yAxe.resize(nCD[1]);
        }
        //Filling the grid
        for (auto particle: particles) {
            int x_cell = (int) (particle->getPosition().getX() / rCut);
            int y_cell = (int) (particle->getPosition().getY() / rCut);
            if (!(x_cell < nCD[0] && y_cell < nCD[1])){
                throw invalid_argument("the particle position surpasses the dimensions of the universe");
            }
            grid[x_cell][y_cell].particles.insert(particle);//2D
            particle->gridPosition.resize(2);
            particle->gridPosition[0] = x_cell;
            particle->gridPosition[1] = y_cell;
        }
    }
    else if (dim == 1){
        //grid
        //Init the grid to the proper size x = nCD[0], y = nCD[1], z = nCD[2]
        this->grid1D.resize(nCD[0]);
        //Filling the grid
        for (auto particle: particles) {
            int x_cell = (int) (particle->getPosition().getX() / rCut);
            if (!(x_cell < nCD[0])){
                throw invalid_argument("the particle position x surpasses the dimension of the universe");
            }

            grid1D[x_cell].particles.insert(particle);//2D
            particle->gridPosition.resize(1);
            particle->gridPosition[0] = x_cell;
        }
    }
    else if (dim == 3){
        //grid
        //Init the grid to the proper size x = nCD[0], y = nCD[1], z = nCD[2]
        this->grid3D.resize(nCD[0]);
        for (auto &yAxe: grid3D) {
            yAxe.resize(nCD[1]);
            for (auto &zAxe : yAxe){
                zAxe.resize(nCD[2]);
            }
        }
        //Filling the grid
        for (auto particle: particles) {
            int x_cell = (int) (particle->getPosition().getX() / rCut);
            int y_cell = (int) (particle->getPosition().getY() / rCut);
            int z_cell = (int) (particle->getPosition().getZ() / rCut);
            if (!(x_cell < nCD[0] && y_cell < nCD[1] && z_cell < nCD[2])){
                throw invalid_argument("the particle position surpasses the dimensions of the universe");
            }
            grid3D[x_cell][y_cell][z_cell].particles.insert(particle);//2D
            particle->gridPosition.resize(3);
            particle->gridPosition[0] = x_cell;
            particle->gridPosition[1] = y_cell;
            particle->gridPosition[2] = z_cell;
        }
    }
    else {
        throw invalid_argument("The dimension entered for the universe is not valid");
    }
};


std::vector<std::vector<Cell>> Universe::getGrid()
{
    return grid;
}



ostream &operator<<(std::ostream &os, const Universe &universe)
{

    if (universe.dim == 2) {
        for (int x = 0; x < universe.nCD[0]; x++) {
            for (int y = 0; y < universe.nCD[1]; y++) {
                os << "Cell: " << x << " , " << y << " particles : ";
                for (auto particle: universe.grid[x][y].getParticles()) {
                    os << *particle << " ; ";
                }
                os << endl;
            }
        }
    }
    else if (universe.dim == 3) {
        for (int x = 0; x < universe.nCD[0]; x++) {
            for (int y = 0; y < universe.nCD[1]; y++) {
                for (int z = 0; z < universe.nCD[2]; z ++) {
                    os << "Cell: " << x << " , " << y << "," << z << " particles : ";
                    for (auto particle: universe.grid3D[x][y][z].getParticles()) {
                        os << *particle << " ; ";
                    }
                    os << endl;
                }
            }
        }
    }
    else if (universe.dim == 1) {
        for (int x = 0; x < universe.nCD[0]; x++) {
            os << "Cell: " << x <<  " particles : ";
            for (auto particle: universe.grid1D[x].getParticles()) {
                os << *particle << " ; ";
            }
            os << endl;
        }
    }

    return os;
}


/**
 * Calculates the potential interaction forces between a particle and its neighboring particles 2D version
 * @param universe
 */
void interactionForcesPotentiel(Universe & universe, bool ljReflexion, double G)
{
    if (universe.boundCond != 0 && ljReflexion){
        throw invalid_argument("The Lennard-Jones boundary condition for reflexion cannot be applied unless boundCond"
                               "for the universe is 0 (Absorption)");
    }
    for (int x = 0; x < universe.nCD[0]; x++){
        for (int y = 0; y < universe.nCD[1]; y++){
            //We are in the cell (x,y)
            //We loop over the neighboring cells, if they exist we consider the particles inside
            for (auto particleI : universe.grid[x][y].getParticles()){
                particleI->force = Vecteur(0.0, 0.0, 0.0);
                if (ljReflexion) {
                    //Detecting the closest boundary and choosing the direction
                    int dirX;double rX;
                    if (particleI->position.getX() < universe.lD[0] - particleI->position.getX()){
                        //Left boundary the velocity is negative, the magnitude should then be positive
                        rX = particleI->position.getX();
                        dirX = 1;
                    }
                    else {
                        rX = universe.lD[0] - particleI->position.getX();
                        dirX  = - 1;
                    }

                    int dirY;double rY;
                    if (particleI->position.getY() < universe.lD[1] - particleI->position.getY()){
                        //Left boundary the velocity is negative, the magnitude should then be positive
                        rY = particleI->position.getY();
                        dirY = 1;
                    }
                    else {
                        rY = universe.lD[1] - particleI->position.getY();
                        dirY  = - 1;
                    }
                    

                    double rCutRef = pow(2, 1.0 / 6.0) * universe.sigma;
                    if (rX < rCutRef && rX != 0){
                        //The direction depends on the speed direction, it is the opposite
                        double flJMagnitude = -24 * universe.epsilon  * (1/(2*rX)) * fastPow(universe.sigma / (2*rX), 6) *
                            (1 - 2 * fastPow(universe.sigma / (2*rX), 6));
                        particleI->force.setX(particleI->force.getX() + dirX * flJMagnitude);
                    }
                    if (rY < rCutRef && rY != 0){
                        double flJMagnitude = -24 * universe.epsilon  * (1/(2*rY)) * fastPow(universe.sigma / (2*rY), 6) *
                                              (1 - 2 * fastPow(universe.sigma / (2*rY), 6));
                        particleI->force.setY(particleI->force.getY() + dirY * flJMagnitude);
                    }
                }

                if(G > 0){
                    particleI->force.setY(particleI->force.getY() - particleI->masse * G);
                }

                //Interactions with same cell particles
                for (auto particleJ : universe.grid[x][y].getParticles()){
                    if (!(particleI == particleJ)){
                        Vecteur rIJVect = particleJ->position - particleI->position;
                        double rIJ = rIJVect.norm();
                        if (rIJ < universe.rCut && rIJ > 0){
                            double fIJMagnitude = 24 * universe.epsilon * (1/(rIJ*rIJ)) *
                                                  fastPow((universe.sigma /rIJ), 6) * (1- 2 * fastPow((universe.sigma/rIJ), 6));
                            particleI->force = particleI->force + rIJVect * fIJMagnitude;
                        }
                    }
                }

                //Interactions with the neighbors
                for (int dx = -1; dx <= 1; dx ++){
                    for (int dy = -1; dy <= 1; dy ++){
                        if (dx == 0 && dy == 0) continue;
                        if (x + dx < universe.nCD[0] && x + dx >= 0 && y + dy < universe.nCD[1] && y + dy >= 0){
                            for (auto particleJ : universe.grid[x + dx][y + dy].getParticles()){
                                Vecteur rIJVect = particleJ->position - particleI->position;
                                double rIJ = rIJVect.norm();
                                if (rIJ < universe.rCut && rIJ > 0){
                                    double fIJMagnitude = 24 * universe.epsilon * (1/(rIJ*rIJ)) *
                                                          fastPow((universe.sigma /rIJ), 6) * (1- 2 * fastPow((universe.sigma/rIJ), 6));
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

/**
 * Calculates the potential interaction forces between a particle and its neighboring particles in 3D
 * @param universe
 */
void interactionForcesPotentiel3D(Universe & universe, bool ljReflexion, double G)
{
    if (universe.boundCond != 0 && ljReflexion){
        throw invalid_argument("The Lennard-Jones boundary condition for reflexion cannot be applied unless boundCond"
                               "for the universe is 0 (Absorption)");
    }
    for (int x = 0; x < universe.nCD[0]; x++){
        for (int y = 0; y < universe.nCD[1]; y++){
            for (int z = 0; z < universe.nCD[2]; z++){
                //We are in the cell (x,y,z)
                for (auto particleI : universe.grid3D[x][y][z].getParticles()){
                    particleI->force = Vecteur(0.0, 0.0, 0.0);
                    if (ljReflexion) {
                        int dirX;double rX;
                        if (particleI->position.getX() < universe.lD[0] - particleI->position.getX()){
                            rX = particleI->position.getX();
                            dirX = 1;
                        }
                        else {
                            rX = universe.lD[0] - particleI->position.getX();
                            dirX  = - 1;
                        }

                        int dirY;double rY;
                        if (particleI->position.getY() < universe.lD[1] - particleI->position.getY()){
                            rY = particleI->position.getY();
                            dirY = 1;
                        }
                        else {
                            rY = universe.lD[1] - particleI->position.getY();
                            dirY  = - 1;
                        }

                        int dirZ; double rZ;
                        if (particleI->position.getZ() < universe.lD[2] - particleI->position.getZ()){
                            rZ = particleI->position.getZ();
                            dirZ = 1;
                        }
                        else {
                            rZ = universe.lD[2] - particleI->position.getZ();
                            dirZ  = - 1;
                        }

                        double rCutRef = pow(2, 1.0 / 6.0) * universe.sigma;
                        if (rX < rCutRef && rX > 0){
                            double flJMagnitude = -24 * universe.epsilon  * (1/(2*rX)) * fastPow(universe.sigma / (2*rX), 6) *
                                                  (1 - 2 * fastPow(universe.sigma / (2*rX), 6));
                            particleI->force.setX(particleI->force.getX() + dirX * flJMagnitude);
                        }
                        if (rY < rCutRef && rY > 0){
                            double flJMagnitude = -24 * universe.epsilon  * (1/(2*rY)) * fastPow(universe.sigma / (2*rY), 6) *
                                                  (1 - 2 * fastPow(universe.sigma / (2*rY), 6));
                            particleI->force.setY(particleI->force.getY() + dirY * flJMagnitude);
                        }
                        if (rZ < rCutRef && rZ > 0){
                            double flJMagnitude = -24 * universe.epsilon  * (1/(2*rZ)) * fastPow(universe.sigma / (2*rZ), 6) *
                                                  (1 - 2 * fastPow(universe.sigma / (2*rZ), 6));
                            particleI->force.setZ(particleI->force.getZ() + dirZ * flJMagnitude);
                        }
                    }

                    if(G > 0){
                        particleI->force.setZ(particleI->force.getZ() - particleI->masse * G);
                    }

                    for (auto particleJ : universe.grid3D[x][y][z].getParticles()){
                        if (!(particleI == particleJ)){
                            Vecteur rIJVect = particleJ->position - particleI->position;
                            double rIJ = rIJVect.norm();
                            if (rIJ < universe.rCut && rIJ > 0){
                                double fIJMagnitude = 24 * universe.epsilon * (1/(rIJ*rIJ)) *
                                                      fastPow((universe.sigma /rIJ), 6) * (1- 2 * fastPow((universe.sigma/rIJ), 6));
                                particleI->force = particleI->force + rIJVect * fIJMagnitude;
                            }
                        }
                    }

                    for (int dx = -1; dx <= 1; dx ++){
                        for (int dy = -1; dy <= 1; dy ++){
                            for (int dz = -1; dz <= 1; dz ++){
                                if (dx == 0 && dy == 0 && dz == 0) continue;
                                if (x + dx < universe.nCD[0] && x + dx >= 0 && y + dy < universe.nCD[1] && y + dy >= 0 && z + dz < universe.nCD[2] && z + dz >= 0){
                                    for (auto particleJ : universe.grid3D[x + dx][y + dy][z + dz].getParticles()){
                                        Vecteur rIJVect = particleJ->position - particleI->position;
                                        double rIJ = rIJVect.norm();
                                        if (rIJ < universe.rCut && rIJ > 0){
                                            double fIJMagnitude = 24 * universe.epsilon * (1/(rIJ*rIJ)) *
                                                                  fastPow((universe.sigma /rIJ), 6) * (1- 2 * fastPow((universe.sigma/rIJ), 6));
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
}

/**
 * Calculates the potential interaction forces between a particle and its neighboring particles in 1D
 * @param universe
 */
void interactionForcesPotentiel1D(Universe & universe, bool ljReflexion, double G)
{
    if (universe.boundCond != 0 && ljReflexion){
        throw invalid_argument("The Lennard-Jones boundary condition for reflexion cannot be applied unless boundCond"
                               "for the universe is 0 (Absorption)");
    }
    for (int x = 0; x < universe.nCD[0]; x++){
        //We are in the cell (x)
        for (auto particleI : universe.grid1D[x].getParticles()){
            particleI->force = Vecteur(0.0, 0.0, 0.0);
            if (ljReflexion) {
                int dirX; double rX;
                if (particleI->position.getX() < universe.lD[0] - particleI->position.getX()){
                    rX = particleI->position.getX();
                    dirX = 1;
                }
                else {
                    rX = universe.lD[0] - particleI->position.getX();
                    dirX  = - 1;
                }

                double rCutRef = pow(2, 1.0 / 6.0) * universe.sigma;
                if (rX < rCutRef && rX != 0){
                    double flJMagnitude = -24 * universe.epsilon  * (1/(2*rX)) * fastPow(universe.sigma / (2*rX), 6) *
                                          (1 - 2 * fastPow(universe.sigma / (2*rX), 6));
                    particleI->force.setX(particleI->force.getX() + dirX * flJMagnitude);
                }
            }

            //G doesn't do anything

            for (auto particleJ : universe.grid1D[x].getParticles()){
                if (!(particleI == particleJ)){
                    Vecteur rIJVect = particleJ->position - particleI->position;
                    double rIJ = rIJVect.norm();
                    if (rIJ < universe.rCut && rIJ > 0){
                        double fIJMagnitude = 24 * universe.epsilon * (1/(rIJ*rIJ)) *
                                              fastPow((universe.sigma /rIJ), 6) * (1- 2 * fastPow((universe.sigma/rIJ), 6));
                        particleI->force = particleI->force + rIJVect * fIJMagnitude;
                    }
                }
            }

            for (int dx = -1; dx <= 1; dx ++){
                if (dx == 0) continue;
                if (x + dx < universe.nCD[0] && x + dx >= 0){
                    for (auto particleJ : universe.grid1D[x + dx].getParticles()){
                        Vecteur rIJVect = particleJ->position - particleI->position;
                        double rIJ = rIJVect.norm();
                        if (rIJ < universe.rCut && rIJ > 0){
                            double fIJMagnitude = 24 * universe.epsilon * (1/(rIJ*rIJ)) *
                                                  fastPow((universe.sigma /rIJ), 6) * (1- 2 * fastPow((universe.sigma/rIJ), 6));
                            particleI->force = particleI->force + rIJVect * fIJMagnitude;
                        }
                    }
                }
            }

        }
    }
}



/**
 * Called after change of position of particles, it removes a particle from the old cell and puts it in the new cell 2D
 * @param universe
 */
void updateGrid(Universe &universe, Particle * particle)
{
    int x_cell = (int)(particle->getPosition().getX() / universe.rCut);
    int y_cell = (int)(particle->getPosition().getY() / universe.rCut);
    //2D
    universe.grid[particle->gridPosition[0]][particle->gridPosition[1]].removeParticle(particle);
    if (x_cell >= 0 && x_cell < universe.nCD[0] && y_cell >= 0 && y_cell < universe.nCD[1]) {
        universe.grid[x_cell][y_cell].particles.insert(particle);//2D
        particle->gridPosition[0] = x_cell;
        particle->gridPosition[1] = y_cell;
    }
}
/**
 * Called after change of position of particles, it removes a particle from the old cell and puts it in the new cell 3D
 * @param universe
 */
void updateGrid3D(Universe &universe, Particle * particle)
{
    int x_cell = (int)(particle->getPosition().getX() / universe.rCut);
    int y_cell = (int)(particle->getPosition().getY() / universe.rCut);
    int z_cell = (int)(particle->getPosition().getZ() / universe.rCut);
    universe.grid3D[particle->gridPosition[0]][particle->gridPosition[1]][particle->gridPosition[2]].removeParticle(particle);
    if (x_cell >= 0 && x_cell < universe.nCD[0] && y_cell >= 0 && y_cell < universe.nCD[1] && z_cell >= 0 && z_cell < universe.nCD[2]) {
        universe.grid3D[x_cell][y_cell][z_cell].particles.insert(particle);
        particle->gridPosition[0] = x_cell;
        particle->gridPosition[1] = y_cell;
        particle->gridPosition[2] = z_cell;
    }
}

/**
 * Called after change of position of particles, it removes a particle from the old cell and puts it in the new cell 1D
 * @param universe
 */
void updateGrid1D(Universe &universe, Particle * particle)
{
    int x_cell = (int)(particle->getPosition().getX() / universe.rCut);
    universe.grid1D[particle->gridPosition[0]].removeParticle(particle);
    if (x_cell >= 0 && x_cell < universe.nCD[0]) {
        universe.grid1D[x_cell].particles.insert(particle);
        particle->gridPosition[0] = x_cell;
    }
}




void printVtk(vector<Particle *> particleList, ostream & outputStream){
    outputStream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n"
                    "  <UnstructuredGrid>\n"
                    "    <Piece NumberOfPoints=\""<< particleList.size()<< "\" NumberOfCells=\"0\">" << endl;
    outputStream << "<Points>\n"
                    "        <DataArray name=\"Position\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;
    for (auto particle : particleList) outputStream << particle->position << " ";

    outputStream << "\n        </DataArray>\n"
                    "      </Points>" << endl;

    outputStream << "<PointData Vectors=\"vector\">\n"
                    "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (auto particle : particleList) outputStream << particle->vitesse << " ";
    outputStream << "        </DataArray>\n"
                    "        <DataArray type=\"Float32\" Name=\"Masse\" format=\"ascii\">\n";
    for (auto particle : particleList) outputStream << particle->masse << " ";
    outputStream  <<"\n        </DataArray>\n"
                    "      </PointData>\n"
                    "      <Cells>\n"
                    "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n"
                    "        </DataArray>\n"
                    "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n"
                    "        </DataArray>\n"
                    "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n"
                    "        </DataArray>\n"
                    "      </Cells>\n"
                    "    </Piece>\n"
                    "  </UnstructuredGrid>\n"
                    "</VTKFile>" << endl;
}



double kineticEnergy(Universe &universe){
    double e = 0;
    for (auto particle : universe.getParticles()){
        e += particle->masse * particle->vitesse.norm()*particle->vitesse.norm();
    }
    return e/2;
}




/**
 * Uses the stromerVerlet integration method to calculate the new position for the system of particles
 * in the universe. It considers the Lennard-Jones potential interaction between particles. It adapts to the boundary
 * condition of the universe. And if demanded by the user it scales the velocity to limit the kinetic energy in
 * the system
 * @param universe
 * @param tEnd : The duration of the simulation
 * @param deltaT : intergration step
 * @param visual : if true the output of the universe is printed as .vtu file to be visualized
 * @param path : The path for the visualization, only makes sense when visual is true
 * @param ljReflexion : Does reflexion on the boundary using the lennard jones force, it shouldn't be passed as true
 * if the boundary condition of the universe (universe.boundCond) is either -1 or 1
 * @param G : The uniform gravitational potential, if G = 0 the potential is not considered
 * @param eCD : The value used for scaling the kinetic energy for the universe, if eCD = 0 there will be no scaling
 */
void stromerVerletPotential(Universe &universe, double tEnd, double deltaT,
                            bool visual, string path, bool ljReflexion, double G, double eCD)
{
    if (universe.boundCond != 0 && ljReflexion){
        throw invalid_argument("The Lennard-Jones boundary condition for reflexion cannot be applied unless boundCond"
                               "for the universe is 0 (Absorption)");
    }


    // Calculate the initial forces
    vector<Particle *> particleList = universe.particles;
    vector<Vecteur> fOld;
    //Show init state
    ofstream os;
    os.open(path + "sim0.vtu");
    printVtk(universe.particles, os);
    os.close();
    //
    interactionForcesPotentiel(universe, ljReflexion, G);
    for (auto particleI : particleList)
    {
        fOld.push_back(particleI->force);
    }

    double t = 0.0;
    int iter = 0;

    int fps = 30;
    int visStepsPerFrame = round(1.0 / (fps * deltaT));
    int frameCounter = 1;

    while (t < tEnd)
    {
        //Testing if to visualize
        if (visual && iter % visStepsPerFrame == 0){
            ofstream os;
            os.open(path + "sim" + to_string(frameCounter) + ".vtu");
            printVtk(universe.particles, os);
            os.close();
            frameCounter++;
        }

        //Testing if to scale speed in the simulation
        if (eCD > 0 && iter >0 && iter % 1000 == 0){
            double beta = sqrt(eCD/ kineticEnergy(universe));
            for (auto particle: particleList){
                particle->vitesse = particle->vitesse * beta;
            }
        }

        //Main part
        iter ++;
        t += deltaT;
        int index = 0;
        if (universe.boundCond == 0 || universe.boundCond == -1) {
            //absorption
            for (auto particleI: particleList) {
                particleI->position =
                    particleI->position + (particleI->vitesse +
                                           particleI->force * (0.5 / particleI->masse) * deltaT) *
                                          deltaT;
                updateGrid(universe, particleI);
                fOld[index] = particleI->force;
                index++;
            }
        }
        else if (universe.boundCond == 1) {
            //periodic
            for (auto particleI: particleList) {
                //portal_x = new_x % x_end
                particleI->position =
                    (particleI->position + (particleI->vitesse +
                                            particleI->force * (0.5 / particleI->masse) * deltaT) *
                                           deltaT);
                
                double newX = particleI->position.getX();
                double newY = particleI->position.getY();
                if (newX >= universe.lD[0]){
                    newX = fmod(newX, universe.lD[0]);
                }
                else if (newX < 0){
                    newX = 100 - abs(fmod(newX, universe.lD[0]));
                }
                
                if (newY >= universe.lD[1]){
                    newY = fmod(newY, universe.lD[1]);
                }
                else if (newY < 0){
                    newY = 100 - abs(fmod(newY, universe.lD[1]));
                }

                particleI->position = Vecteur(newX, newY, 0.0);

                updateGrid(universe, particleI);

                fOld[index] = particleI->force;
                index++;
            }
        }

        interactionForcesPotentiel(universe, ljReflexion, G);

        index = 0; // will be reused after this
        for (auto particleI : particleList)
        {
            particleI->vitesse = particleI->vitesse + (particleI->force + fOld[index]) *
                                                      deltaT * (0.5 / particleI->masse);
            if (universe.boundCond == -1){
                //Reflexion
                if (particleI->position.getX() >= universe.lD[0] || particleI->position.getX() <= 0){
                    particleI->vitesse.setX( -particleI->vitesse.getX());
                }
                if (particleI->position.getY() >= universe.lD[1] || particleI->position.getY() <= 0){
                    particleI->vitesse.setY( -particleI->vitesse.getY());
                }
            }
            index++;
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

const vector<Particle *> &Universe::getParticles() const {
    return particles;
}

void Universe::setGrid1D(const vector<Cell> &grid1D) {
    Universe::grid1D = grid1D;
}

void Universe::setGrid3D(const vector<std::vector<std::vector<Cell>>> &grid3D) {
    Universe::grid3D = grid3D;
}



/**
 * Uses the stromerVerlet integration method to calculate the new position for the system of particles
 * in the universe. It considers the Lennard-Jones potential interaction between particles. It adapts to the boundary
 * condition of the universe. And if demanded by the user it scales the velocity to limit the kinetic energy in
 * the system
 * @param universe
 * @param tEnd : The duration of the simulation
 * @param deltaT : intergration step
 * @param visual : if true the output of the universe is printed as .vtu file to be visualized
 * @param path : The path for the visualization, only makes sense when visual is true
 * @param ljReflexion : Does reflexion on the boundary using the lennard jones force, it shouldn't be passed as true
 * if the boundary condition of the universe (universe.boundCond) is either -1 or 1
 * @param G : The uniform gravitational potential, if G = 0 the potential is not considered
 * @param eCD : The value used for scaling the kinetic energy for the universe, if eCD = 0 there will be no scaling
 */
void stromerVerletPotential3D(Universe &universe, double tEnd, double deltaT,
                            bool visual, string path, bool ljReflexion, double G, double eCD)
{
    if (universe.boundCond != 0 && ljReflexion){
        throw invalid_argument("The Lennard-Jones boundary condition for reflexion cannot be applied unless boundCond"
                               "for the universe is 0 (Absorption)");
    }

    // Calculate the initial forces
    vector<Particle *> particleList = universe.particles;
    vector<Vecteur> fOld;
    //Show init state
    ofstream os;
    os.open(path + "sim0.vtu");
    printVtk(universe.particles, os);
    os.close();
    //
    interactionForcesPotentiel3D(universe, ljReflexion, G);
    for (auto particleI : particleList)
    {
        fOld.push_back(particleI->force);
    }

    double t = 0.0;
    int iter = 0;

    int fps = 30;
    int visStepsPerFrame = round(1.0 / (fps * deltaT));
    int frameCounter = 1;

    while (t < tEnd)
    {
        //Testing if to visualize
        if (visual && iter % visStepsPerFrame == 0){
            ofstream os;
            os.open(path + "sim" + to_string(frameCounter) + ".vtu");
            printVtk(universe.particles, os);
            os.close();
            frameCounter++;
        }

        //Testing if to scale speed in the simulation
        if (eCD > 0 && iter >0 && iter % 1000 == 0){
            double beta = sqrt(eCD/ kineticEnergy(universe));
            for (auto particle: particleList){
                particle->vitesse = particle->vitesse * beta;
            }
        }

        //Main part
        iter ++;
        t += deltaT;
        int index = 0;
        if (universe.boundCond == 0 || universe.boundCond == -1) {
            //absorption
            for (auto particleI: particleList) {
                particleI->position =
                    particleI->position + (particleI->vitesse +
                                           particleI->force * (0.5 / particleI->masse) * deltaT) *
                                          deltaT;
                updateGrid3D(universe, particleI);
                fOld[index] = particleI->force;
                index++;
            }
        }
        else if (universe.boundCond == 1) {
            //periodic
            for (auto particleI: particleList) {
                //portal_x = new_x % x_end
                particleI->position =
                    (particleI->position + (particleI->vitesse +
                                            particleI->force * (0.5 / particleI->masse) * deltaT) *
                                           deltaT);

                double newX = particleI->position.getX();
                double newY = particleI->position.getY();
                double newZ = particleI->position.getZ();
                if (newX >= universe.lD[0]){
                    newX = fmod(newX, universe.lD[0]);
                }
                else if (newX < 0){
                    newX = 100 - abs(fmod(newX, universe.lD[0]));
                }

                if (newY >= universe.lD[1]){
                    newY = fmod(newY, universe.lD[1]);
                }
                else if (newY < 0){
                    newY = 100 - abs(fmod(newY, universe.lD[1]));
                }

                if (newZ >= universe.lD[2]){
                    newZ = fmod(newZ, universe.lD[2]);
                }
                else if (newZ < 0){
                    newZ = 100 - abs(fmod(newZ, universe.lD[2]));
                }
                particleI->position = Vecteur(newX, newY, newZ);

                updateGrid3D(universe, particleI);

                fOld[index] = particleI->force;
                index++;
            }
        }

        interactionForcesPotentiel3D(universe, ljReflexion, G);

        index = 0; // will be reused after this
        for (auto particleI : particleList)
        {
            particleI->vitesse = particleI->vitesse + (particleI->force + fOld[index]) *
                                                      deltaT * (0.5 / particleI->masse);
            if (universe.boundCond == -1){
                //Reflexion
                if (particleI->position.getX() >= universe.lD[0] || particleI->position.getX() <= 0){
                    particleI->vitesse.setX( -particleI->vitesse.getX());
                }
                if (particleI->position.getY() >= universe.lD[1] || particleI->position.getY() <= 0){
                    particleI->vitesse.setY( -particleI->vitesse.getY());
                }
                if (particleI->position.getZ() >= universe.lD[2] || particleI->position.getZ() <= 0){
                    particleI->vitesse.setZ( -particleI->vitesse.getZ());
                }
            }
            index++;
        }
    }
}

void stromerVerletPotential1D(Universe &universe, double tEnd, double deltaT,
                            bool visual, string path, bool ljReflexion, double G, double eCD)
{
    if (universe.boundCond != 0 && ljReflexion){
        throw invalid_argument("The Lennard-Jones boundary condition for reflexion cannot be applied unless boundCond"
                               "for the universe is 0 (Absorption)");
    }


    // Calculate the initial forces
    vector<Particle *> particleList = universe.particles;
    vector<Vecteur> fOld;
    //Show init state
    ofstream os;
    os.open(path + "sim0.vtu");
    printVtk(universe.particles, os);
    os.close();
    //
    interactionForcesPotentiel(universe, ljReflexion, G);
    for (auto particleI : particleList)
    {
        fOld.push_back(particleI->force);
    }

    double t = 0.0;
    int iter = 0;

    int fps = 30;
    int visStepsPerFrame = round(1.0 / (fps * deltaT));
    int frameCounter = 1;

    while (t < tEnd)
    {
        //Testing if to visualize
        if (visual && iter % visStepsPerFrame == 0){
            ofstream os;
            os.open(path + "sim" + to_string(frameCounter) + ".vtu");
            printVtk(universe.particles, os);
            os.close();
            frameCounter++;
        }

        //Testing if to scale speed in the simulation
        if (eCD > 0 && iter >0 && iter % 1000 == 0){
            double beta = sqrt(eCD/ kineticEnergy(universe));
            for (auto particle: particleList){
                particle->vitesse = particle->vitesse * beta;
            }
        }

        //Main part
        iter ++;
        t += deltaT;
        int index = 0;
        if (universe.boundCond == 0 || universe.boundCond == -1) {
            //absorption
            for (auto particleI: particleList) {
                particleI->position =
                    particleI->position + (particleI->vitesse +
                                           particleI->force * (0.5 / particleI->masse) * deltaT) *
                                          deltaT;
                updateGrid1D(universe, particleI);
                fOld[index] = particleI->force;
                index++;
            }
        }
        else if (universe.boundCond == 1) {
            //periodic
            for (auto particleI: particleList) {
                //portal_x = new_x % x_end
                particleI->position =
                    (particleI->position + (particleI->vitesse +
                                            particleI->force * (0.5 / particleI->masse) * deltaT) *
                                           deltaT);

                double newX = particleI->position.getX();
                if (newX >= universe.lD[0]){
                    newX = fmod(newX, universe.lD[0]);
                }
                else if (newX < 0){
                    newX = 100 - abs(fmod(newX, universe.lD[0]));
                }

                particleI->position = Vecteur(newX, 0.0, 0.0);

                updateGrid1D(universe, particleI);

                fOld[index] = particleI->force;
                index++;
            }
        }

        interactionForcesPotentiel1D(universe, ljReflexion, G);

        index = 0; // will be reused after this
        for (auto particleI : particleList)
        {
            particleI->vitesse = particleI->vitesse + (particleI->force + fOld[index]) *
                                                      deltaT * (0.5 / particleI->masse);
            if (universe.boundCond == -1){
                //Reflexion
                if (particleI->position.getX() >= universe.lD[0] || particleI->position.getX() <= 0){
                    particleI->vitesse.setX( -particleI->vitesse.getX());
                }
            }
            index++;
        }
    }
}







/**
 * It checks the dimension of the universe to call the proper integration function
 * Uses the stromerVerlet integration method to calculate the new positions for the system of particles
 * in the universe. It considers the Lennard-Jones potential interaction between particles. It adapts to the boundary
 * condition of the universe. And if demanded by the user it scales the velocity to limit the kinetic energy in
 * the system
 * @param universe
 * @param tEnd : The duration of the simulation
 * @param deltaT : intergration step
 * @param visual : if true the output of the universe is printed as .vtu file to be visualized
 * @param path : The path for the visualization, only makes sense when visual is true
 * @param ljReflexion : Does reflexion on the boundary using the lennard jones force, it shouldn't be passed as true
 * if the boundary condition of the universe (universe.boundCond) is either -1 or 1
 * @param G : The uniform gravitational potential, if G = 0 the potential is not considered
 * @param eCD : The value used for scaling the kinetic energy for the universe, if eCD = 0 there will be no scaling
 */
void startSimulation(Universe &universe, double tEnd, double deltaT,
                     bool visual, string path, bool ljReflexion, double G, double eCD){
    if (universe.boundCond != 0 && ljReflexion){
        throw invalid_argument("The Lennard-Jones boundary condition for reflexion cannot be applied unless boundCond"
                               "for the universe is 0 (Absorption)");
    }

    if (universe.dim == 2){
        stromerVerletPotential(universe, tEnd, deltaT, visual, path, ljReflexion, G, eCD);
    }
    else if (universe.dim == 1){
        stromerVerletPotential1D(universe, tEnd, deltaT, visual, path, ljReflexion, G, eCD);
    }
    else if (universe.dim == 3){
        stromerVerletPotential3D(universe, tEnd, deltaT, visual, path, ljReflexion, G, eCD);
    }
}



