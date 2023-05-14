#include "vecteur.hxx"


bool Vecteur::operator<(const Vecteur &v2) const{
        return (*this).norm() < v2.norm();
}

bool Vecteur::operator==(const Vecteur &v2){
        return x == v2.x && y  == v2.y && z == v2.z;
}

ostream& operator<<(ostream& o, const Vecteur &v){
        return o << v.x << "    " << v.y << "    " << v.z;
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

