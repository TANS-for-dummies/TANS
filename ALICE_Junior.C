//Librerie custom
#include "MyRandom.h"
#include "Rivelatore.h"
#include "Particella.h"
#include "Vertice.h" //Carica anche Punto.h

//Librerie
#include "Riostream.h"
#include "TMath.h"

void MonteCarlo(int gen = 1, unsigned int seed = 125) {
    
    //Costanti
    double pi_greco = TMath::Pi();
    
    //Settaggi
    const char* input_file = "kinem.root";
    int N = 30;
    
    //Generatore di numeri random
	MyRandom *ptr = new MyRandom(input_file,seed);
	delete gRandom;
	gRandom = ptr;
	
    //Creazione del funtore per scegliere la molteplicita'
    int (MyRandom::*rndm_molt) (int); 
    if(gen == 1) {rndm_molt = &MyRandom::RndMolt;}
    else if (gen == 2) {rndm_molt = &MyRandom::RndMolt_unif;}
    else if (gen == 3) {rndm_molt = &MyRandom::RndMolt_fissa;}
    else {cout << "Scelta non valida. Impostato il settaggio di base: estrazione dall'istogramma" << endl;
        rndm_molt = &MyRandom::RndMolt;}

	//Rivelatori
	Rivelatore Beam_Pipe(3,0.08,4);
	Rivelatore Layer1(4,0.02,14);
	Rivelatore Layer2(7,0.02,14);
    
    if (ptr->GetFlag()) {
    	std::cout << "File con le distribuzioni non trovato" << endl; 
    	return;
    }

    //Iniziamo a generare il vertice, ci servono 3 coordinate e la molteplicità
    //Unità di misura della lunghezza = cm
    Vertice Ver(ptr->Gaus(0.,0.01),ptr->Gaus(0.,0.01),ptr->Gaus(0.,5.3),(ptr->*rndm_molt)(N));
    
    //Generiamo i prodotti nel vertice
    for(int i=0;i<Ver.GetN();i++) {
    	Ver.AddPart(ptr->RndTheta(),ptr->Rndm()*2.*pi_greco);
    	cout << "Particella Numero " << i+1 << " : ( " << Ver.GetPart(i).GetTheta() << " , " <<  Ver.GetPart(i).GetPhi() << " )" << endl; 
    }
   
}
