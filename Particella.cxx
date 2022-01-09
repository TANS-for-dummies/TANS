#include "Particella.h"
#include "TMath.h"

//classe importata da root
ClassImp(Particella)

//Costruttore di default
Particella::Particella() : dmCoord1(0.), dmCoord2(0.){}

//Costruttore standard
Particella::Particella(double C1, double C2): dmCoord1(C1), dmCoord2(C2), TObject(){}

//Costruttore d copia
Particella::Particella(const Particella& source) : TObject(source){
    dmCoord1 = source.dmCoord1;
    dmCoord2 = source.dmCoord2;
}

//Distruttore standard
Particella::~Particella(){}

//Operatore copia
Particella& Particella::operator=(const Particella& source){
    if(this == &source) return *this;
    this->~Particella();
    new(this) Particella(source);
    return *this;
}

