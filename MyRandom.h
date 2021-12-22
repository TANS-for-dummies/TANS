#ifndef MYRANDOM_H
#define MYRANDOM_H

#include "TRandom3.h"
#include "TH1D.h"
#include "TFile.h"

class MyRandom : public TRandom3 {

    public:

        MyRandom(); //Costruttore di default
        MyRandom(const char* input_file, double seed); //Costruttore standard
        MyRandom(const MyRandom& source); //Copy-Constructor
        virtual ~MyRandom(); //Distruttore
        MyRandom& operator=(const MyRandom& source); //Operatore =
        
        //Estraiamo dagli istogrammi valori casuali
        double RndTheta1();
        int RndMolt() {return (int) dmMolt->GetRandom();};
        
        static bool GetFlag() {return sFlag;}
        
    private:
    	static bool sFlag;
    	
    	TH1D* dmEta;
    	TH1D* dmMolt;
	
    ClassDef(MyRandom,1)
};

#endif
