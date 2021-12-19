#include "Rivelatore.h"

ClassImp(Rivelatore)

//Costruttore di default
Rivelatore::Rivelatore(): dmR(0.), dmS(0.), dmZ(1), TObject() {}

//Costruttore standard
Rivelatore::Rivelatore(double r, double s, int z): dmR(r), dmS(s), dmZ(z), TObject() {}

//Copy
Rivelatore::Rivelatore(const Rivelatore& source) : dmR(source.dmR), dmS(source.dmS), dmZ(source.dmZ), TObject(source) {}

//Distruttore
Rivelatore::~Rivelatore {}

//Operatore copia
Rivelatore& Rivelatore::operator=(const Rivelatore& source){
    if(this == &source) return *this;
    this->~Rivelatore();
    new(this) Rivelatore(source);
    return *this;
}
