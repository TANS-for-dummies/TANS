//Librerie custom
#include "MyRandom.h"
#include "Rivelatore.h"
#include "Particella.h"
#include "Segnale.h"
#include "Punto.h"

//Librerie
#include "Riostream.h"
#include "TMath.h"
#include "TStopwatch.h" //Monitora il tempo di CPU
#include "TSystem.h" //Monitora RAM
#include "TTree.h" //Output
#include "TBranch.h"
#include "TLeaf.h"
#include "TClonesArray.h"
#include <fstream>



void MonteCarlo(int N_esp = 1000000, int gen = 1, bool scat = 1, unsigned int seed = 125) {
	
    //Costanti
    double pi_greco = TMath::Pi();
    
    //Settaggi input e output
    int N_false_hit = 0;
    const char* input_file = "kinem.root";
    int N = 30;
    const char* output_file = "MonteCarlo.root";
    //double Theta_Multi = (1/TMath::Sqrt2())*0.001; //rad
    double Theta_Multi_Be = 13.6*4*TMath::Sqrt(0.08*1.85/65.19)*(1+0.038*TMath::Log(0.08*1.85/65.19))/(3.*pow(10,13));
    double Theta_Multi_Si = 13.6*14*TMath::Sqrt(0.02*2.33/22.0)*(1+0.038*TMath::Log(0.02*2.33/22.0))/(3.*pow(10,13));

    //Avviamo il timer	
    TStopwatch timer;
    timer.Start();
    
    //Apro il file di log
    std::ofstream ofs ("simulation_log.txt", std::ofstream::out);

	
    //Per monitorare la RAM
    ProcInfo_t* proc = new ProcInfo_t();
    ofs << "RAM utilizzata: " << gSystem->GetProcInfo(proc) << std::endl;

    //Generatore di numeri random
    MyRandom *ptr = new MyRandom(input_file,seed);
    delete gRandom;
    gRandom = ptr;

    if (ptr->GetFlag()) {
    	ofs << "File con le distribuzioni non trovato" << std::endl; 
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
    else {ofs << "Scelta non valida. Impostato il settaggio di base: estrazione dall'istogramma" << std::endl;
        rndm_molt = &MyRandom::RndMolt;
        dim = 36;
        }



    //Creazione del funtore per MultiScattering
    Particella (Rivelatore::*rndm_scatt) (Particella*, MyRandom*);
    if (scat){rndm_scatt = &Rivelatore::MultiScattering;}
    else {rndm_scatt = &Rivelatore::ZeroScattering;}

    
    //Apertura file di output, e creazione di un TTree
    TFile Ofile(output_file, "RECREATE");
    TTree *tree = new TTree("Tree","TTree con 3 branches"); //Vertice, layer1 e layer2

    TClonesArray *riv_1 = new TClonesArray("Segnale",dim);//Hit del rivelatore 1
    TClonesArray &hit1 = *riv_1;
    TClonesArray *riv_2 = new TClonesArray("Segnale",dim);//Hit del rivelatore 2
    TClonesArray &hit2 = *riv_2;
    


	//Rivelatori
	Rivelatore Beam_Pipe(3, 0.08,52, Theta_Multi_Be); //H=52 per contenere tutte le particelle generate con l'accettanza data
	Rivelatore Layer1(4, 0.02, 27, Theta_Multi_Si);
	Rivelatore Layer2(7, 0.02, 27, Theta_Multi_Si);
    

    // Definiamo una struct 
    typedef struct{
        double x, y, z;
        int molt;
    } Vertice;
    
    static Vertice inizio;
    
    tree->Branch("VertMult", &inizio.x, "x/D:y:z:molt/I"); 
    tree->Branch("Hit1", &riv_1);
    tree->Branch("Hit2", &riv_2);

    tree->SetAutoSave(0); //Rimuove backup cycle
    

    //Creiamo una particella (fuori dal for così viene creata una sola volta), sarà descritta da 2 angoli
    Particella* part = new Particella();

    //Creiamo un hit temporaneo
    Punto* hit = new Punto();
    
    
    for(int k=0; k<N_esp; k++){
        //Iniziamo a generare il vertice, ci servono 3 coordinate e la molteplicità
        //Unità di misura della lunghezza = cm
        inizio.x = ptr->Gaus(0.,0.01);
        inizio.y = ptr->Gaus(0.,0.01);
        inizio.z = ptr->Gaus(0.,5.3);
        inizio.molt = (ptr->*rndm_molt)(N);
        
	if (k/1000 == 0){ofs << "(" << inizio.x << ", "<<inizio.y << ", "<<inizio.z << ") e molteplicità " << inizio.molt << std::endl;};
	
        int pos1 = 0;
        int pos2 = 0;

        for(int i=0; i<inizio.molt; i++) {
            //Generiamo i prodotti nel vertice
            part->SetTheta(ptr->RndTheta());
            part->SetPhi(ptr->Rndm()*2.*pi_greco);
            if (k/1000 == 0){ofs << "Particella Numero " << i+1 << " : ( " << part->GetTheta() << " , " <<  part->GetPhi() << " )" << std::endl;}; 
            
            //Trasporto e multiscattering particella per particella

            

            //BEAM PIPE
            *hit = Beam_Pipe.Hit(Punto(inizio.x, inizio.y, inizio.z), part);
            *part = (Beam_Pipe.*rndm_scatt)(part, ptr);
            
            //LAYER 1
            *hit = Layer1.Hit(*hit, part);
            if (TMath::Abs(hit -> GetZ())>((Layer1.GetH())/2.)){}
            else{

                //Immagazziniamo lo smearing (coordinale cilindriche)
                new(hit1[pos1]) Segnale(Layer1.Smearing(hit, ptr, i+1));

                //Multiscattering
                *part = (Layer1.*rndm_scatt)(part, ptr);

                //LAYER 2
                *hit = Layer2.Hit(*hit, part);
                if(TMath::Abs(hit -> GetZ())>((Layer2.GetH())/2.)){}
                else{
                    //Immagazziniamo lo smearing
                    new(hit2[pos2]) Segnale(Layer2.Smearing(hit, ptr, i+1));
                    pos2++;
                }
                pos1++;
            }

        }

        //Generazione dei false hit (uniformi in Z e phi)
        for(int i=0; i<N_false_hit; i++) {
            new(hit1[pos1+i]) Segnale(-(Layer1.GetH())/2.+(ptr->Rndm())*(Layer1.GetH()),ptr->Rndm()*2*pi_greco, -(i+1));
            new(hit2[pos2+i]) Segnale(-(Layer2.GetH())/2.+(ptr->Rndm())*(Layer2.GetH()),ptr->Rndm()*2*pi_greco, -(i+1));
        }

        // Debug
        printf("Entries nel TClonesArray: %d\n",riv_1->GetEntries());
        for (int j=0; j<hit1.GetEntries(); j++){
        Particella *tst = (Particella*)hit1[j];
        // Particella *tst=(Particella*)riv_1->At(j);
        if (k/1000 == 0){ofs <<"Particella "<<j+1<<") Z , phi = "<<tst->GetTheta()<<"; "<<tst->GetPhi()<<std::endl;};
        //delete tst;
        }
        printf("Entries nel TClonesArray: %d\n",riv_2->GetEntries());
        for (int j=0; j<hit2.GetEntries(); j++){
        Particella *tst = (Particella*)hit2[j];
        //Particella *tst=(Particella*)riv_2->At(j);
        if (k/1000 == 0){ofs <<"Particella "<<j+1<<") Z , phi = "<<tst->GetTheta()<<"; "<<tst->GetPhi()<<std::endl;};
        //delete tst;
        }
        // fine del debug
        

        tree->Fill();

        riv_1->Clear();
        riv_2->Clear();
    }

    // Save all objects in this file
    Ofile.Write();
    
    // Close the file. 
    Ofile.Close();
    
    timer.Stop();
    timer.Print();

    ofs << "RAM utilizzata: " << gSystem->GetProcInfo(proc) << std::endl;
	
    ofs.close();
    delete part;
    delete hit;

}
