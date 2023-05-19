#include "cell.hxx"

void Cell::removeParticle(Particle * particle){
    particles.erase(particle);
}

const unordered_set<Particle *> &Cell::getParticles() const {
    return particles;
}
