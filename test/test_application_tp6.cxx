#include <iostream>
#include <random>
#include <list>

#include "universe.hxx"
#include "particule.hxx"




int main() {


    double lD[2] = {250, 180};
    double epsilon = 1;
    double sigma = 1;
    double rCut = 2.5*sigma;
    double d_t = 0.00005;

    vector<Particle * > particleList;

    //Cercle
    int centerX = (int)round(lD[0] / 2); // Circle center X
    int centerY = 101; // Circle center Y
    int radius = 15; // Circle radius (12% of the caracteristic distance
    int N1 = 395; // Number of particles

    double spacingCercle = sqrt(M_PI * radius * radius / N1); // Average spacing between the particles

    // Calculate the number of rings and number of particles per ring
    int numRings = radius / spacingCercle;
    int numParticlesPerRing = N1 / numRings;

    for (int i = 0; i < numRings; i++) {
        double ringRadius = i * spacingCercle;
        int numParticlesThisRing = numParticlesPerRing;
        // Special case: the last ring takes all remaining particles
        if (i == numRings - 1)
            numParticlesThisRing = N1 - i * numParticlesPerRing;
        double angleIncrement = 2.0 * M_PI / numParticlesThisRing;
        for (int j = 0; j < numParticlesThisRing; j++) {
            // Calculate angle
            double angle = j * angleIncrement;
            // Calculate particle's position
            double x = centerX + ringRadius * cos(angle);
            double y = centerY + ringRadius * sin(angle);
            string type = "nil";
            Particle *particleToInsert = new Particle(Vecteur(x, y, 0.0),
                                                      Vecteur(0.0, 10.0, 0.0), 1.0, i * numParticlesPerRing + j, Vecteur(0.0, 0.0, 0.0), type);
            particleList.push_back(particleToInsert);
        }
    }

    //Rectangle
    int N2 = 17227; // Number of particles
    int length = 250; // Length of the rectangle
    int height = 72; // Height of the rectangle

    double area = (double)(length * height); // Area of the rectangle
    double spacing = sqrt(area / N2); // Spacing between the particles

    int numParticlesX = (int)(length / spacing); // Number of particles along X
    int numParticlesY = (int)(height / spacing); // Number of particles along Y

    for (int i = 0; i < numParticlesX; i++) {
        for (int j = 0; j < numParticlesY; j++) {
            // Calculate particle's position
            double x = i * spacing;
            double y = j * spacing;

            string type = "nil";
            Particle *particleToInsert = new Particle(Vecteur(x, y, 0.0),
                                                      Vecteur(0.0, 10.0, 0.0), 1.0, i * numParticlesY + j, Vecteur(0.0, 0.0, 0.0), type);
            particleList.push_back(particleToInsert);
        }
    }

    // Handle remaining particles
    int remainingParticles = N2 - numParticlesX * numParticlesY;
    for (int i = 0; i < remainingParticles; i++) {
        double x = (i % numParticlesX) * spacing;
        double y = (i / numParticlesX) * spacing;

        string type = "nil";
        Particle *particleToInsert = new Particle(Vecteur(x, y, 0.0),
                                                  Vecteur(0.0, 10.0, 0.0), 1.0, numParticlesX * numParticlesY + i, Vecteur(0.0, 0.0, 0.0), type);
        particleList.push_back(particleToInsert);
    }




    Universe *univ = new Universe(particleList, 2, rCut, lD, epsilon, sigma);

    ofstream init;
    init.open("../../demo/initTP6.vtu");
    printVtk(univ->getParticles(), init);

//    stromerVerletPotential(*univ, 29.5, d_t, init, "../../demo/");

//    printVtk(univ->getParticles(), init);
    init.close();
//    end.close();

    return 0;
};