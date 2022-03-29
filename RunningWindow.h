#ifndef RUNNINGWINDOW_H
#define RUNNINGWINDOW_H

#include "TObject.h"

class RunningWindow : public TObject {

	public:		
		
		//Costruttori e associati vari da inserire
		RunningWindow(){
            dmSize = 0.;
            dmStep = 1.;
        };
		RunningWindow(double size, double step){
            dmSize = size;
            dmStep = step;
        };
		RunningWindow(const RunningWindow& source);
		virtual ~RunningWindow(){};
		RunningWindow& operator=(const RunningWindow& source);		
		
		/*
		I getter, fanno il return di un data member sono dei const 		
		siccome non devono poter alterare il punto
		*/
		double GetSize() const {return dmSize;};
        double GetStep() const {return dmStep;};

		
		//I setter, impostano il valore di uno specifico data member
		void SetSize(double size) {dmSize=size;};
        void SetStep(double step) {dmStep=step;};

        //Funzione di analisi
        double running_window(vector<double> vec,bool &stato_rec);

        //media degli elementi del vector
        double media(vector<double> V,int j,double limite);

        static void ResetRaddoppio(){raddoppio = 0;};
	
	private:
	
        static bool raddoppio;

		//Data member
        double dmSize;
        double dmStep;

	
	//Definisce questa come versione 1  di RunningWindow per ROOT
	ClassDef(RunningWindow,1)
	
};

#endif