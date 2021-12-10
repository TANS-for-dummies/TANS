#include "Particella.h"

ClassImp(Particella)

//Costruttore di default
Particella::Particella() : dmTheta(0.), dmPhi(0.), dmP(0.){}

//Costruttore standard
Particella::Particella(double Theta, double Phi, double P): dmTheta(Theta), dmPhi(Phi), dmP(P), TObject(){}

//Costruttore d copia
Particella::Particella(const Particella& source) : TObject(source){
    dmTheta = source.dmTheta;
    dmPhi = source.dmPhi;
    dmP = source.dmP;
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