#ifndef PARTICLE_H
#define PARTICLE_H

class Universe;

#include <string>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>
#include <set>
#include "vecteur.hxx"


using namespace std;


class Particle {
private:
        Vecteur position;
        vector<int> gridPosition;
        Vecteur vitesse;
        double  masse;

private:
        int identifiant;
        Vecteur force;
        string type;
public:
        Particle(){};
        Particle(Vecteur position, Vecteur vitesse, double masse, int identifiant, Vecteur force, string type):
        position(position), vitesse(vitesse),
         masse(masse), identifiant(identifiant), force(force), type(type){
        }


        Vecteur getPosition() const;
        /**members**/
        bool operator<(const Particle&) const;
        bool operator==(const Particle&);

        /**
         * Return the Lennard Jones potential between "this" particle and another particle
         * p is defined as a pointer because it is normally allocated dynamically
         * @param p
         * @param epsilon
         * @param sigma
         * @return
         */
        double interactionLJ(Particle *p,double epsilon, double sigma);




        /**friends**/
        friend void stromerVerlet(vector<Particle*> &particleList, double tEnd, double deltaT, ofstream &outputStream);
        friend void calculateForces(vector<Particle*> &particleList);
        friend void calculateForcesSlow(vector<Particle*> &particleList);
        friend ostream& operator<<(ostream &o, const Particle &);
        friend void printVtk(vector<Particle *> particleList, ostream & outputStream);

        friend void interactionForcesPotentiel(Universe & universe, bool ljRelexion, double G);
        friend void interactionForcesPotentiel3D(Universe & universe, bool ljRelexion, double G);
        friend void interactionForcesPotentiel1D(Universe & universe, bool ljRelexion, double G);
        friend void stromerVerletPotential(Universe & universe, double tEnd, double deltaT,bool visual,  string path,
                                           bool ljReflexion, double G, double eCD);
        friend void stromerVerletPotential3D(Universe & universe, double tEnd, double deltaT,bool visual,  string path,
                                           bool ljReflexion, double G, double eCD);
        friend void stromerVerletPotential1D(Universe & universe, double tEnd, double deltaT,bool visual,  string path,
                                           bool ljReflexion, double G, double eCD);

        friend double kineticEnergy(Universe &universe);
        friend void updateGrid(Universe &universe, Particle * particle);
        friend void updateGrid1D(Universe &universe, Particle * particle);
        friend void updateGrid3D(Universe &universe, Particle * particle);



        friend class Universe;

        void setPosition(const Vecteur &position);
};



/**
 * Simulate the movement of the particles using their inter-forces
 * @param particleList
 * @param tEnd
 * @param deltaT
 * @param outputStream
 */
void stromerVerlet(vector<Particle*> &particleList, double tEnd, double deltaT, ofstream &outputStream, bool ljReflexion);




/**
 * Calculate the interaction forces in particleList
 * @param particleList : reference to a vector of pointers to particle instances
 */
void calculateForces(vector<Particle *> &particleList);






list<Particle> constructParticleList(int n, bool print);

set<Particle> constructParticleSet(int n, bool print);

vector<Particle> constructParticleVector(int n, bool print);

deque<Particle> constructParticleDeque(int n, bool print);

#endif