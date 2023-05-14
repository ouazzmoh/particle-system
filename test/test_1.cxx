#include <iostream>

#include "vecteur.hxx"


int main(){
        Vecteur v1;
        Vecteur v2;
        Vecteur v3 = Vecteur(3.0, 2.0, 1.0);
        Vecteur v4 = Vecteur(2,1,1);

        // test <<
        std::cout<<v1<<std::endl;
        std::cout<<v2<<std::endl;

        // test norm
        if (v2 < v3){
                cout << "The norm of" << v3 << " is :" << v3.norm() <<"\n";
        }

        // test operator +
        std::cout<< "\ntest operator + : \n"<< v3 << " + " << v4 << " = " <<  v3 + v4<<std::endl;

        // test operator +=
        Vecteur v5 ;
        v5 = v4;
        v5 += v3;

        std::cout<< "\ntest operator += : \n"<< "v5 = " <<  v4 << "\n" << "v5 += v3  ==> v5  = " << v5 << std::endl;

        // test operator *
        std::cout<< "\ntest operator * : \n"<< v3 << " * " << v4 << " = " <<  v3 * v4<<std::endl;

        // test operator *=
        Vecteur v6 ;
        v6 = v4;
        v6 *= v3;
        std::cout<< "\ntest operator *= : \n"<< "v6 = " <<  v4 << "\n" << "v6 *= v3  ==> v6  = " << v6 << std::endl;

        // test operator -
        std::cout<< "\ntest operator - : \n"<< v3 << " - " << v4 << " = " <<  v3 - v4<<std::endl;

        // test operator -=
        Vecteur v7 ;
        v7 = v4;
        v7 -= v3;
        std::cout<< "\ntest operator -= : \n"<< "v7 = " <<  v4 << "\n" << "v7 -= v3  ==> v7  = " << v7 << std::endl;

        // test operator * double
        std::cout<< "\ntest operator * : \n"<< v3 << " * 4" << " = " <<  v3 * 4<<std::endl;

        // test operator *=
        Vecteur v8 ;
        v8 = v3;
        v8 *= 4;
        std::cout<< "\ntest operator *= for double : \n"<< "v8 = " <<  v3 << "\n" << "v8 *= 4  ==> v8  = " << v8 << std::endl;

        // test  <
        bool b = v3 < v4;
        std::cout<< "\nnorme de v3  : " << v3.norm() << "\nnorme de v4 : " << v4.norm() << "\nv3  < v4 = "  <<   b <<std::endl;
        // test  ==
        bool b1 = v3 == v4;
        std::cout<< "\nv3 =  " << v3<< "\nv4 = " << v4 << "\nb1 = v3  == v4 ==> b1 = "  <<   b1  <<std::endl;
        bool b2 = v3 == v3;
        std::cout<< "\nv3 =  " << v3 << "\nb2  = v3  == v3  ==> b2 = "  <<   b2  <<std::endl;




    return EXIT_SUCCESS;
}