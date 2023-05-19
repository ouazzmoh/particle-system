#ifndef VECTEUR_H
#define VECTEUR_H

#include <cmath>
#include <ostream>

using namespace std;


/**
 * Faster function to calculate power of integers
 * @param base
 * @param exponent
 * @return
 */
double fastPow(double base, int exponent);


class Vecteur {
        double x, y, z;
public:

        Vecteur(double x = 0.0, double y = 0.0, double z = 0.0):x(x), y(y), z(z){}


        double norm() const{
                return sqrt(x*x + y*y + z*z);
        };
        /**
         * Returns the distance to another point in polar cordinates ("this" is the reference)
         * We assume the axis are x and y

         * @param v
         * @return
         */
        double distanceToVect(Vecteur v) const{
                return sqrt((x - v.x)*(x - v.x) + (y - v.y)*(y - v.y) + (z - v.z)*(z - v.z));
        }



        /**Comparaison functions**/
        bool operator<(const Vecteur &v2) const;
        bool operator==(const Vecteur &v2);




        /**Operations*/
        Vecteur operator+(const Vecteur &v2);
        Vecteur& operator+=(const Vecteur &v2);
        Vecteur operator*(const Vecteur &v2);
        Vecteur& operator*=(const Vecteur &v2);
        Vecteur operator-(const Vecteur &v2);
        Vecteur& operator-=(const Vecteur &v2);
        Vecteur& operator%(const Vecteur &v2);

        Vecteur operator*(double f);
        Vecteur& operator*=(double f);


        /**Showing vectors**/
        friend ostream& operator<<(ostream&, const Vecteur &v);

        double getX() const;

        double getY() const;

        double getZ() const;

        void setX(double x);

        void setY(double y);

        void setZ(double z);


};

ostream& operator<<(ostream&, const Vecteur &v);




#endif