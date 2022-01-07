//File di implementazione della classe Punto

#include "Punto.h"

#include "Riostream.h"
#include "TMath.h"

ClassImp(Punto)

//Costruttore di default
Punto::Punto() : dmX(0), dmY(0), dmZ(0) {
	//Per ora va bene quello standard
}

//Costruttore standard
Punto::Punto(double x, double y, double z) : dmX(x), dmY(y), dmZ(z), TObject() {}

//Copy-Constructor
Punto::Punto(const Punto& source) : TObject(source) {
	//Non allocando memoria possiamo usarne uno semplice
	dmX=source.dmX;
	dmY=source.dmY;
	dmZ=source.dmZ;
}

//Destructor
Punto::~Punto() {
	//Per ora va bene quello standard
}

//Operatore =
Punto& Punto::operator=(const Punto& source) {
	//Utilizzo quello che sfrutta il copy constructor per non ripetere tutto
	if(this == &source) return *this;
	this->~Punto();
	new(this) Punto(source);
	return *this;
}

//I Getter e i setter li ho già definiti nel .h

//GetRadiusXY
double Punto::GetRadiusXY() const {
	return TMath::Sqrt(dmX*dmX+dmY*dmY);
}
