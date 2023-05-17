#include "cell.hxx"

void Cell::removeParticle(Particle * particle){
    particles.erase(std::remove(particles.begin(), particles.end(), particle), particles.end());
}

const vector<Particle *> &Cell::getParticles() const {
    return particles;
}
