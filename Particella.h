/*La classe particella descrive una particella tramite 2 coordinate eventualmente 2 angoli o 1 angolo+1 lunghezza*/

#ifndef PARTICELLA_H
#define PARTICELLA_H

#include "TObject.h"
#include "Punto.h"


class Particella: public TObject{

    public:
        //Costruttore standard
        Particella();

        //Costruttore di default
        Particella(double C1, double C2);

        //Copy
        Particella(const Particella& source);

        //Distruttore: virtual poichè viene utilizzato il distruttore di TObject (e non alloca memoria)
        virtual ~Particella();

        //Overloading operatore =
        Particella& operator=(const Particella& source);

        //GETTER
        double GetCoord1() const {return dmCoord1;};
        double GetCoord2() const {return dmCoord2;};

        //SETTER (bau)
        void SetCoord1(double C1){ dmCoord1 = C1;};
        void SetCoord2(double C2){ dmCoord2 = C2;};

        



    private:
        //Data members
        double dmCoord1;
        double dmCoord2;
        


ClassDef(Particella,1)
};
#endif
