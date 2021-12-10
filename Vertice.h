/*La classe serve a descrivere l'oggetto VERTICE come: 
punto di interazione, direzione delle particelle prodotte e molteplicità.*/


#ifndef VERTICE_H
#define VERTICE_H

#include "Punto.h"
#include "Particella.h"
#include "vector"

class Vertice: public Punto{

    public:
        //Costruttore standard
        Vertice();

        //Costruttore di default
        Vertice(double x, double y, double z, int N);

        //Copy
        Vertice(const Vertice& source);

        //Distruttore
        ~Vertice();

        //Overloading operatore =
        Vertice& operator=(const Vertice& source);

        //GETTER
        int GetN() const {return dmN;};
        Particella GetPart(int i) const {return dmDir.At(i);}; 

        //SETTER
        void SetN(double N) {dmN = N;};

        void AddPart(double Theta, double Phi, double p);


    private:
        //Data members
        int dmN;
        std::vector<Particella> dmDir;


ClassDef(Vertice,1)
};
#endif
