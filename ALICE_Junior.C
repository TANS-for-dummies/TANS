//Librerie custom
#include "MyRandom.h"
#include "Rivelatore.h"
#include "Particella.h"
#include "Punto.h" 

//Librerie
#include "Riostream.h"
#include "TMath.h"
#include <TStopwatch.h> //Monitora il tempo di CPU
#include "TSystem.h" //Monitora RAM
#include "TTree.h" //Output
#include "TBranch.h"
#include "TClonesArray.h"

void MonteCarlo(int gen = 1, bool scat = 1, unsigned int seed = 125) {
	
    //Costanti
    double pi_greco = TMath::Pi();
    
    //Settaggi input e output
    const char* input_file = "kinem.root";
    int N = 30;
    const char* output_file = "MonteCarlo.root";
    //double Theta_Multi = (1/TMath::Sqrt2())*0.001; //rad
    double Theta_Multi_Be = 13.6*4*TMath::Sqrt(0.08*1.85/65.19)*(1+0.038*TMath::Log(0.08*1.85/65.19))/(3.*pow(10,13));
    double Theta_Multi_Si = 13.6*14*TMath::Sqrt(0.02*2.33/22.0)*(1+0.038*TMath::Log(0.02*2.33/22.0))/(3.*pow(10,13));

    //Avviamo il timer	
    TStopwatch timer;
    timer.Start();
	
    ProcInfo_t* proc = new ProcInfo_t();
    cout << "RAM utilizzata: " << gSystem->GetProcInfo(proc) << endl;

    //Generatore di numeri random
    MyRandom *ptr = new MyRandom(input_file,seed);
    delete gRandom;
    gRandom = ptr;

    if (ptr->GetFlag()) {
    	std::cout << "File con le distribuzioni non trovato" << endl; 
    	return;
    }


    //Creazione del funtore per scegliere la molteplicita'
    int dim = 0;
    int (MyRandom::*rndm_molt) (int); 
    if(gen == 1) {
        rndm_molt = &MyRandom::RndMolt;
        dim = 36; //68.27% di 53 (massimo valore della molteplicità)
        }
    else if (gen == 2) {
        rndm_molt = &MyRandom::RndMolt_unif;
        dim = N/2 +1;
        }
    else if (gen == 3) {
        rndm_molt = &MyRandom::RndMolt_fissa;
        dim = N;
        }
    else {cout << "Scelta non valida. Impostato il settaggio di base: estrazione dall'istogramma" << endl;
        rndm_molt = &MyRandom::RndMolt;
        dim = 36;
        }



    //Creazione del funtore per MultiScattering
    Particella (Rivelatore::*rndm_scatt) (Particella*, MyRandom*);
    if (scat){rndm_scatt = &Rivelatore::MultiScattering;}
    else {rndm_scatt = &Rivelatore::ZeroScattering;}

    /*
    //Apertura filedi output, e creazione di un TTree
    TFile Ofile(output_file, "RECREATE");
    TTree *tree = new TTree("Tree","TTree con 3 branches"); //Vertice, layer1 e layer2

    TClonesArray *riv_1 = new TClonesArray("Coord_cil",dim);//Hit del rivelatore 1
    TClonesArray &hit1 = *riv_1;
    TClonesArray *riv_2 = new TClonesArray("Coord_cil",dim);//Hit del rivelatore 2
    TClonesArray &hit2 = *riv_2;
    */


	//Rivelatori
	Rivelatore Beam_Pipe(3, 0.08,52, Theta_Multi_Be); //H=52 per contenere tutte le particelle generate con l'accettanza data
	Rivelatore Layer1(4, 0.02, 27, Theta_Multi_Si);
	Rivelatore Layer2(7, 0.02, 27, Theta_Multi_Si);
    

    // Definiamo una struct 
    typedef struct{
    Punto P;
    int molt;} Vertice;
    
    static Vertice inizio;

    static Coord_cil segnale;
    /*
    tree->Branch("VertMult", &inizio.P, "P.dmX/D:P.dmY:P.dmZ:molt/I"); //DA GUARDARE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    tree->Branch("Hit1", &riv_1);
    tree->Branch("Hit2", &riv_2);
    */

    //Creiamo una particella (fuori dal for così viene creata una sola volta)
    Particella* part = new Particella();

    //Creiamo un hit temporaneo
    Punto* hit = new Punto();

    //Iniziamo a generare il vertice, ci servono 3 coordinate e la molteplicità
    //Unità di misura della lunghezza = cm
    inizio.P = Punto(ptr->Gaus(0.,0.01),ptr->Gaus(0.,0.01),ptr->Gaus(0.,5.3));
    inizio.molt = (ptr->*rndm_molt)(N);
    
    cout<< "(" << inizio.P.GetX()<<", "<<inizio.P.GetY()<<", "<<inizio.P.GetZ()<<")"<<endl;


    
    for(int i=0; i<inizio.molt; i++) {
        //Generiamo i prodotti nel vertice
        part->SetTheta(ptr->RndTheta());
        part->SetPhi(ptr->Rndm()*2.*pi_greco);
    	cout << "Particella Numero " << i+1 << " : ( " << part->GetTheta() << " , " <<  part->GetPhi() << " )" << endl; 
        //Trasporto e multiscattering particella per particella

        //BEAM PIPE
        *hit = Beam_Pipe.Hit(inizio.P, part);
        *part = (Beam_Pipe.*rndm_scatt)(part, ptr);
        
        //LAYER 1
        *hit = Layer1.Hit(*hit, part);
        if (TMath::Abs(hit -> GetZ())>Layer1.GetH()){}
        else{

            //Immagazziniamo lo smearing
            //segnale = Layer1.Smearing(hit, ptr);
            //new(hit1[i]) Coord_cil(Layer1.Smearing(hit, ptr));

            //Multiscattering
            *part = (Layer1.*rndm_scatt)(part, ptr);

            //LAYER 2
            *hit = Layer2.Hit(*hit, part);
            if(TMath::Abs(hit -> GetZ())>Layer2.GetH()){}
            else{
                //Immagazziniamo lo smearing
                //segnale = Layer2.Smearing(hit, ptr);
                //new(hit2[i]) Coord_cil(Layer2.Smearing(hit, ptr));
            }

        }

    }
    /*
    tree->Fill();
    riv_1->Clear();
    riv_2->Clear();

    // Save all objects in this file
    Ofile.Write();
    
    // Close the file. 
    Ofile.Close();
    */

    
	
    timer.Stop();
    timer.Print();

    cout << "RAM utilizzata: " << gSystem->GetProcInfo(proc) << endl;
	
    delete part;
    delete hit;

}
