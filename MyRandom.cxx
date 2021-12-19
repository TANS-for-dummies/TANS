#include "MyRandom.h"
#include "Riostream.h"
#include "TMath.h"

ClassImp(MyRandom)

//Costruttore di default
MyRandom::MyRandom() : TRandom3() {}

//Costruttore standard
MyRandom::MyRandom(double seed) : TRandom3(seed) {}

//Copy-Constructor
MyRandom::MyRandom(const MyRandom& source) : TRandom3(source) {}

//Destructor
MyRandom::~MyRandom() {}

//Operatore =
MyRandom& MyRandom::operator=(const MyRandom& source) {
    if(this == &source) return *this;
    this->~MyRandom();
    new(this) MyRandom(source);
    return *this;
}
