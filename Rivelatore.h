/*
La classe rivelatore si occupa di fornire i due layer sensibili ma anche la
beam pipe
*/

#ifndef RIVELATORE_H
#define RIVELATORE_H

#include "TObject.h"

class Rivelatore : public TObject {

	public:
		//Costruttore standard
		Rivelatore();
		
		//Costruttore di default
		Rivelatore(double r, double s, int z);
		
		//Copy
		Rivelatore(const Rivelatore& source);
		
		//Distruttore
		virtual ~Rivelatore();
		
		//Overloading dell'operatore =
		Rivelatore& operator=(const Rivelatore& source);
		
		//GETTER
		double GetR() const {return dmR;};
		double GetS() const {return dmS;};
		int GetZ() const {return dmZ;};
		
		//SETTER
		void SetR(double r) {dmR=r;};
		void SetS(double s) {dmS=s;};
		void SetZ(int z) {dmZ=z;};
	
	private:
		double dmR; //Raggio
		double dmS; //Spessore
		int dmZ;    //Z del materiale

ClassDef(Rivelatore,1)
};
#endif
