#include "vecteur.hxx"
#include <cmath>


bool Vecteur::operator<(const Vecteur &v2) const{
        return (*this).norm() < v2.norm();
}

bool Vecteur::operator==(const Vecteur &v2){
        return x == v2.x && y  == v2.y && z == v2.z;
}

ostream& operator<<(ostream& o, const Vecteur &v){
        return o << v.x << " " << v.y << " " << v.z;
}



Vecteur Vecteur::operator+(const Vecteur &v2){
        return Vecteur(x + v2.x, y + v2.y, z + v2.z);
}
Vecteur & Vecteur::operator+=(const Vecteur &v2){
        x += v2.x;
        y += v2.y;
        z += v2.z;
        return *this;
}
Vecteur Vecteur::operator*(const Vecteur &v2){
        return Vecteur(x * v2.x, y * v2.y, z * v2.z);
}
Vecteur & Vecteur::operator*=(const Vecteur &v2){
        x *= v2.x;
        y *= v2.y;
        z *= v2.z;
        return *this;
}
Vecteur Vecteur::operator-(const Vecteur &v2){
        return Vecteur(x - v2.x, y - v2.y, z - v2.z);
}


Vecteur & Vecteur::operator-=(const Vecteur &v2){
        x -= v2.x;
        y -= v2.y;
        z -= v2.z;
        return *this;
}


Vecteur & Vecteur::operator%(const Vecteur &v2){
    if (x > 0){
        x = abs(fmod(x, v2.x));
    }
    else {

    }

    x = abs(fmod(x, v2.x));
    y = abs(fmod(y, v2.y));
    z = abs(fmod(z, v2.z));
    return *this;
}


Vecteur Vecteur::operator*(double f){
        return Vecteur(x * f, y * f, z * f);
}
Vecteur & Vecteur::operator*=(double f){
        x *= f;
        y *= f;
        z *= f;
        return *this;
}

double Vecteur::getX() const {
        return x;
}

double Vecteur::getY() const {
        return y;
}

double Vecteur::getZ() const {
        return z;
}



void Vecteur::setX(double x) {
    Vecteur::x = x;
}

void Vecteur::setY(double y) {
    Vecteur::y = y;
}

void Vecteur::setZ(double z) {
    Vecteur::z = z;
}

