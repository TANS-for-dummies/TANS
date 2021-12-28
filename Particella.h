/*La classe particella descrive la direzione di una particella*/

#ifndef PARTICELLA_H
#define PARTICELLA_H

#include "TObject.h"

class Particella: public TObject{

    public:
        //Costruttore standard
        Particella();

        //Costruttore di default
        Particella(double Theta, double Phi);

        //Copy
        Particella(const Particella& source);

        //Distruttore: virtual poichè viene utilizzato il distruttore di TObject (e non alloca memoria)
        virtual ~Particella();

        //Overloading operatore =
        Particella& operator=(const Particella& source);

        //GETTER
        double GetTheta() const {return dmTheta;};
        double GetPhi() const {return dmPhi;};

        //SETTER (bau)
        void SetTheta(double Theta){ dmTheta = Theta;};
        void SetPhi(double Phi){ dmPhi = Phi;};



    private:
        //Data members
        double dmTheta;
        double dmPhi;
        


ClassDef(Particella,1)
};
#endif
