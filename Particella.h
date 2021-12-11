/*La classe particella descrive la direzione ed il modulo del momento*/

#ifndef PARTICELLA_H
#define PARTICELLA_H

#include "TObject.h"

class Particella: public TObject{

    public:
        //Costruttore standard
        Particella();

        //Costruttore di default
        Particella(double Theta, double Phi, double P);

        //Copy
        Particella(const Particella& source);

        //Distruttore
        ~Particella();

        //Overloading operatore =
        Particella& operator=(const Particella& source);

        //GETTER
        double GetTheta() const {return dmTheta;};
        double GetPhi() const {return dmPhi;};
        double GetP() const {return dmP;};

        //SETTER (bau)
        void SetTheta(double Theta){ dmTheta = Theta;};
        void SetPhi(double Phi){ dmPhi = Phi;};
        void SetP(double P){ dmP = P;};



    private:
        //Data members
        double dmTheta;
        double dmPhi;
        double dmP;
        


ClassDef(Particella,1)
};
#endif
