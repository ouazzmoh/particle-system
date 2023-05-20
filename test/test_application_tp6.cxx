#include <iostream>
#include <random>
#include <list>

#include "universe.hxx"
#include "particule.hxx"



using namespace std;
int main() {


    double lD[2] = {408, 180};
    double epsilon = 1;
    double sigma = 1;
    double rCut = 2.5*sigma;
    double d_t = 0.00005;
    double G = 12;
    double eCD = 0.05;

    vector<Particle * > particleList;

    //Circle
    int centerX = (int)round(lD[0] / 2); // Circle center X
    int centerY = 101; // Circle center Y
    int radius = 15;
    int N1 = 395; // Number of particles

    // Minimum distance between particles
    double minDistance = pow(2.0, 1.0/6.0)/ sigma;

    // Spacing between the particles (minimum distance)
    double spacingCercle = minDistance ;

    // Calculate the number of rings
    int numRings = radius / spacingCercle;

    // Initialize a counter for total particles
    int totalParticles = 0;
    int id = 0;

    for (int i = 0; i < numRings; i++) {
        double ringRadius = i * spacingCercle;

        // Calculate number of particles that can fit on this ring
        int numParticlesThisRing = round(2.0 * M_PI * ringRadius / spacingCercle);

        //The last ring takes all remaining particles, but do not exceed N1
        if (i == numRings - 1 || totalParticles + numParticlesThisRing > N1) {
            numParticlesThisRing = N1 - totalParticles;
        }

        double angleIncrement = 2.0 * M_PI / numParticlesThisRing;
        for (int j = 0; j < numParticlesThisRing; j++) {
            // Calculate angle
            double angle = j * angleIncrement;
            // Calculate particle's position
            double x = centerX + ringRadius * cos(angle);
            double y = centerY + ringRadius * sin(angle);
            string type = "nil";
            Particle *particleToInsert = new Particle(Vecteur(x, y, 0.0),
                                                      Vecteur(0.0, -10.0, 0.0), 1.0, i * numParticlesThisRing + j, Vecteur(0.0, 0.0, 0.0), type);
            particleList.push_back(particleToInsert);
            id ++;
        }
        totalParticles += numParticlesThisRing;

        // If we've reached the maximum number of particles, break the loop
        if (totalParticles >= N1) {
            break;
        }
    }


    int N2 = 17227;
    double spacingRectangle = minDistance;
    //Rectangle
    for (int i = 0; i < 364; i++){
        for (int j = 0; j < 47; j ++) {
            string type = "nil";
            Particle *particleToInsert = new Particle(Vecteur(i * spacingRectangle, j * spacingRectangle, 0.0),
                                                      Vecteur(0.0, 0.0, 0.0), 1.0, id, Vecteur(0.0, 0.0, 0.0),
                                                      type);
            particleList.push_back(particleToInsert);
            id++;
        }

    }

    Universe *univ = new Universe(particleList, 2, rCut, lD, epsilon, sigma);
    startSimulation(*univ, 29.5, 0.00005, true, "../../demo/tp6_application/", true, G, eCD);

    return 0;
};