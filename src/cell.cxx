#include "cell.hxx"

using namespace std;

void Cell::removeParticle(Particle * particle){
    particles.erase(particle);
}

const unordered_set<Particle *> &Cell::getParticles() const {
    return particles;
}