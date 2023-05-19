#include <iostream>
#include <random>
#include <list>
#include <chrono>
#include <set>
#include <cmath>
#include <vector>
#include <fstream>

#include "particule.hxx"
#include "vecteur.hxx"
#include "universe.hxx"
#include "cell.hxx"
using namespace std;

/**
 * On compare les points par la distance euclidienne
 * @param p1
 * @param p2
 * @return
 */
bool Particle::operator<(const Particle& p2) const{
        return position < p2.position;
}

/**
 * Deux particules égaux si ils ont le meme identifiant
 * @param p1
 * @param p2
 * @return
 */
bool Particle::operator==(const Particle& p2){
        return identifiant == p2.identifiant;
}

/**
 * Afficher les particules
 * @param o
 * @param p
 * @return
 */
ostream& operator<<(ostream &o, const Particle &p){
    return o << "Id : " << p.identifiant <<  " ---- " << "Position : " << p.position << " ---- "
        "Speed : " << p.vitesse << " ---- "
        "Masse : " << p.masse << " ---- "
        "Force : " << p.force << " ---- "
        "Type : " << p.type;
}




void calculateForces(vector<Particle *> & particleList) {
        //The method below divides time complexity by half
        //Store all the forces when calculating to reuse the calculated values
        Vecteur forces[particleList.size()];
        for (int i = 0; i<particleList.size(); i++){forces[i] = Vecteur(0,0,0);}
        //Iterators for i and j iterations over the particleList
        //We increment the iterator so it points to the first element
        //The position of the iterator should match in current index
        auto itParticleI = particleList.begin();

        for (size_t i = 0; i < particleList.size(); i++) {
                auto itParticleJ = particleList.begin();
                std::advance(itParticleJ, i + 1);

                for (size_t j = i + 1; j < particleList.size(); j++) {
                        double rIJ = (*itParticleI)->getPosition().distanceToVect((*itParticleJ)->getPosition());
                        Vecteur fIJ = ((*itParticleJ)->getPosition() - (*itParticleI)->getPosition()) * (((*itParticleI)->masse * (*itParticleJ)->masse) / pow(rIJ, 3));
                        forces[i] += fIJ;
                        // fJI = -fIJ
                        forces[j] += fIJ * (-1);
                        itParticleJ++;
                }
                (*itParticleI)->force = forces[i];
                itParticleI++;
        }
        
        /*
                vector<Particle *>::iterator itParticleI;
        vector<Particle *>::iterator itParticleJ;
        itParticleI = particleList.begin();

        for (int i = 0; i<particleList.size(); i++){
                itParticleJ = particleList.begin();
                advance(itParticleJ, i+1);
                for (int j = i+1; j < particleList.size(); j++){
                        double rIJ = (*itParticleI)->getPosition().distanceToVect((*itParticleJ)->getPosition());
                        Vecteur fIJ = ((*itParticleJ)->getPosition() - (*itParticleI)->getPosition()) * (((*itParticleI)->masse * (*itParticleJ)->masse) / pow(rIJ, 3));
                        forces[i] += fIJ;
                        //fJI = -fIJ
                        forces[j] += fIJ * (-1); //
                        itParticleJ++;
                }
                (*itParticleI)->force = forces[i];
                itParticleI++;
        }
    }
    */
}

void calculateForcesSlow(vector<Particle *> & particleList) {
        for (auto particleI : particleList){
                Vecteur fI = Vecteur(0.0, 0.0, 0.0);
                for (auto particleJ : particleList){
                        if (!(*particleI == *particleJ)){
                                double rIJ = particleI->getPosition().distanceToVect(particleJ->getPosition());
                                // f12 = ((x2 - x1)i + (y2-y1)j + (z2 - z1)k) m1*m2 / r12^3   (We express the vector r12 in cartesian coords)
                                fI = fI + (particleJ->getPosition() - particleI->getPosition()) * ((particleI->masse * particleJ->masse) / pow(rIJ, 3));
                        }
                }
                particleI->force = fI;
        }
}


/**
 * Simulate the movement of the particles
 * @param particleList
 * @param tEnd
 * @param deltaT
 * @param outputStream
 */
void stromerVerlet(vector<Particle*> &particleList, double tEnd, double deltaT,
                   ofstream &outputStream) {
        //Calculate the initial forces
        vector<Vecteur> fOld;
        calculateForces(particleList);
        for (auto &particleI: particleList) {
                fOld.push_back(particleI->force);
        }
        double t = 0.0;
        while (t < tEnd) {
                t += deltaT;
                int index = 0;
                for (auto &particleI: particleList) {
                        particleI->position =
                                particleI->position + ( particleI->vitesse +
                                particleI->force * (0.5 / particleI->masse) * deltaT) * deltaT;
                        fOld[index] = particleI->force;
                        index++;
                }
                calculateForces(particleList);
                index = 0; //will be reused after this
                for (auto &particleI: particleList) {
                        particleI->vitesse = particleI->vitesse + (particleI->force + fOld[index]) *
                                deltaT * (0.5/particleI->masse);
                        index++;
                }
                for (auto &particleI: particleList) {
                        outputStream << particleI->position << endl;
                }
        }
}



list<Particle> constructParticleList(int n, bool print){
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0.0, 100.0);

    list<Particle> particleList;

    for (int i = 0; i < n; ++i){
        string type = "(nil)";
        Particle particleToInsert = Particle(Vecteur(dist(mt), dist(mt), dist(mt)), Vecteur(0.0, 0.0, 0.0)
                                             , 0.0, i, Vecteur(0.0, 0.0, 0.0), type);
        particleList.insert(particleList.end(), particleToInsert);
    }

    if (print){
        for (auto  &particle : particleList){
            cout << particle << endl;
        }
    }

    return particleList;
}


set<Particle> constructParticleSet(int n, bool print){
        random_device rd;
        mt19937 mt(rd());
        uniform_real_distribution<double> dist(0.0, 100.0);

        set<Particle> particleSet;

        for (int i = 0; i < n; ++i){
                string type = "(nil)";
                Particle particleToInsert = Particle(Vecteur(dist(mt), dist(mt), dist(mt)), Vecteur(0.0, 0.0, 0.0)
                        , 0.0, i, Vecteur(0.0, 0.0, 0.0), type);
                particleSet.insert(particleToInsert);
        }

        if (print){
                for (auto  &particle : particleSet){
                        cout << particle << endl;
                }
        }
        return particleSet;
}

vector<Particle> constructParticleVector(int n, bool print){
        random_device rd;
        mt19937 mt(rd());
        uniform_real_distribution<double> dist(0.0, 100.0);

        vector<Particle> particleVect;

        for (int i = 0; i < n; ++i){
                string type = "(nil)";
                Particle particleToInsert = Particle(Vecteur(dist(mt), dist(mt), dist(mt)), Vecteur(0.0, 0.0, 0.0)
                        , 0.0, i, Vecteur(0.0, 0.0, 0.0), type);
                particleVect.push_back(particleToInsert);
        }

        if (print){
                for (auto  &particle : particleVect){
                        cout << particle << endl;
                }
        }
        return particleVect;
}

deque<Particle> constructParticleDeque(int n, bool print){
        random_device rd;
        mt19937 mt(rd());
        uniform_real_distribution<double> dist(0.0, 100.0);

        deque<Particle> particleDeque;

        for (int i = 0; i < n; ++i){
                string type = "(nil)";
                Particle particleToInsert = Particle(Vecteur(dist(mt), dist(mt), dist(mt)), Vecteur(0.0, 0.0, 0.0)
                        , 0.0, i, Vecteur(0.0, 0.0, 0.0), type);
                particleDeque.push_back(particleToInsert);
        }

        if (print){
                for (auto  &particle : particleDeque){
                        cout << particle << endl;
                }
        }
        return particleDeque;        
}



double Particle::interactionLJ(Particle *p, double epsilon, double sigma){
        double r = position.distanceToVect(p->getPosition());
        assert( r > 0);
        double quotient  = sigma/(position.distanceToVect(p->getPosition()));
        return  4 * epsilon * pow(quotient, 6.0) * ( pow(quotient, 6.0) -  1);
}

Vecteur Particle::getPosition() const {
        return position;
}


vector<long> Particle::calculateGridIndex(double rCut){
        vector<long> result;
        result.push_back((long)(position.getX() / rCut));
        result.push_back((long)(position.getY() / rCut));
        result.push_back((long)(position.getZ() / rCut));
        return result;
}

bool Particle::isCellPositionChanged() const {
    return cellPositionChanged;
}

void Particle::setCellPositionChanged(bool cellPositionChanged) {
    Particle::cellPositionChanged = cellPositionChanged;
}

void Particle::setPosition(const Vecteur &position) {
    Particle::position = Vecteur(position.getX(), position.getY(), position.getZ());
}


void printParticleList(const vector<Particle *> &particleList){
    for (auto  &particle : particleList){
        cout << *particle << endl;
    }
}