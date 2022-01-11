/*
La classe rivelatore si occupa di fornire i due layer sensibili ma anche la
beam pipe
*/

#ifndef RIVELATORE_H
#define RIVELATORE_H

#include "TObject.h"
#include "Particella.h"
#include "Segnale.h"
#include "Punto.h"
#include "MyRandom.h"

class Rivelatore : public TObject {

	public:
		//Costruttore standard
		Rivelatore();
		
		//Costruttore di default
		Rivelatore(double r, double s, double H, double Theta);
		
		//Copy
		Rivelatore(const Rivelatore& source);
		
		//Distruttore
		virtual ~Rivelatore();
		
		//Overloading dell'operatore =
		Rivelatore& operator=(const Rivelatore& source);
		
		//METODI
		Particella ZeroScattering(Particella *part, MyRandom *ptr) {return *part;};
		Particella MultiScattering(Particella *part, MyRandom *ptr);

		Segnale Smearing(Punto *P, MyRandom *ptr, int Num_part);

		Punto Hit(Punto P, Particella *part);


		//GETTER
		double GetR() const {return dmR;};
		double GetS() const {return dmS;};
		double GetH() const {return dmH;};
		double GetTheta() const {return dmTheta;};
		
		//SETTER
		void SetR(double r) {dmR=r;};
		void SetS(double s) {dmS=s;};
		void SetH(double H) {dmH=H;};
		void SetTheta(double Theta) {dmTheta=Theta;};
	
	private:
		double dmR; //Raggio
		double dmS; //Spessore
		double dmH; //Lunghezza 
		double dmTheta;    //Angolo planare medio di multiscattering

ClassDef(Rivelatore,1)
};
#endif
