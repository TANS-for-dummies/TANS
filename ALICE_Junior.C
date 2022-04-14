//Librerie custom
#include "Rivelatore.h" //le altre classi (MyRandom, Particella, Segnale e Punto) sono incluse in Rivelatore

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



void MonteCarlo(int N_esp = 1000000, const char* output_file = "MonteCarlo.root", int gen_molt = 1, int gen_z = 1, bool scat = 1, bool smear = 1, int N_false_hit = 0, unsigned int seed = 125, int verboseEvent = 3154) {
    /*
    N_esp è il numero di esperimenti che si vuole simulare
    output_file è il nome del file di output che verrà generato dal programma, l'estensione deve sempre essere .root
    gen_molt seleziona il tipo di generazione della molteplicità (distribuzione data=1, uniforme=2, costante=3)
    gen_z seleziona il tipo di generazione delle z (distribuzione gaussiana=1, uniforme=2, costante=3)
    scat attiva (1) o disattiva (0) il multiscattering sui rivelatori
    smear attiva (1) o disattiva (0) lo smearing degli hit
    N_false_hit indica il numero di hit falsi che devono essere generati uniformemente sul rivelatore
    */
	
    //Costanti
    double pi_greco = TMath::Pi();
    double Theta_Multi = 0.001/(TMath::Sqrt2()); //rad, deviazione standard della gaussiana che utilizziamo che estrarre il theta' del multiscattering
    
    //Settaggi input e output
    const char* input_file = "kinem.root"; //file con le distribuzioni di molteplicità e pseudorapidità assegnate
  

    //Avviamo il timer	
    TStopwatch timer;
    timer.Start();
    
    //Apro il file di log
    std::ofstream ofs ("simulation_log.txt", std::ofstream::out);

    //Generatore di numeri random
    MyRandom *ptr = new MyRandom(input_file,seed);
    delete gRandom;
    gRandom = ptr;

    //Controlliamo che il file di input con le distribuzioni venga caricato correttamente
    if (ptr->GetFlag()) {
    	ofs << "File con le distribuzioni non trovato" << std::endl; 
    	return;
    }
	
	//Creiamo una variabile per fissare la z nel caso di generazione secondo distribuzione costante
	double z_fissa = 0.;
	
	if(gen_z==3){
		std::cout << "Valore coordinata z del vertice:" << std::endl;    
		std::cin >> z_fissa;
		std::cout << std::endl; 
	}

	
    //Creiamo una struct per indicare come sono state generate la molteplicità e la z del vertice
	typedef struct : TObject{
		int Generazione_z;
        double z_costante;
        int Generazione_molt;
    } Generazione;
    
    Generazione settaggi;
	
	settaggi.Generazione_z = gen_z;
	settaggi.z_costante = z_fissa;
		
	
    //Creazione del funtore per scegliere la molteplicita'
    int dim = 0;
    int N;
    int (MyRandom::*rndm_molt) (int); 
    if(gen_molt == 1) { //distribuzione estratta da grafico fornito
        rndm_molt = &MyRandom::RndMolt;
        dim = 36; //68.27% di 53 (massimo valore della molteplicità)
		settaggi.Generazione_molt = 0;
        }
    else if (gen_molt == 2) { //distribuzione uniforme
	std::cout << "Numero massimo di particelle generabile con distribuzione uniforme:" << std::endl;    
	std::cin >> N;
	std::cout << std::endl;    
        rndm_molt = &MyRandom::RndMolt_unif;
        dim = N/2 +1;
		settaggi.Generazione_molt = -N;
        }
    else if (gen_molt == 3) { //distribuzione fissa
	std::cout << "Numero di particelle da generare:" << std::endl;    
	std::cin >> N;
	std::cout << std::endl; 
        rndm_molt = &MyRandom::RndMolt_fissa;
        dim = N;
		settaggi.Generazione_molt = N;
        }
    else {ofs << "Scelta non valida. Impostato il settaggio di base: estrazione dall'istogramma" << std::endl;
        rndm_molt = &MyRandom::RndMolt;
        dim = 36;
		settaggi.Generazione_molt = 0;  
        }
    


    //Creazione del funtore per MultiScattering e per lo smearing
    Particella (Rivelatore::*rndm_scatt) (Particella*);
    if (scat){rndm_scatt = &Rivelatore::MultiScattering;}
    else {rndm_scatt = &Rivelatore::ZeroScattering;}

    Segnale (Rivelatore::*smearing) (Punto*, int);
    if(smear){smearing = &Rivelatore::Smearing;}
    else {smearing = &Rivelatore::NoSmearing;}

    
    //Apertura file di output, e creazione di un TTree
    TFile Ofile(output_file, "RECREATE");
	settaggi.Write("Generazione");
	
    TTree *tree = new TTree("Tree","TTree con 3 branches"); //Vertice, layer1 e layer2

    TClonesArray *riv_1 = new TClonesArray("Segnale",dim);//Hit del rivelatore 1
    TClonesArray &hit1 = *riv_1;
    TClonesArray *riv_2 = new TClonesArray("Segnale",dim);//Hit del rivelatore 2
    TClonesArray &hit2 = *riv_2;
    
    //Rivelatori
    Rivelatore Beam_Pipe(3, 0.08, 52, Theta_Multi, 0., 0.); //H=52 per contenere tutte le particelle generate con l'accettanza data
    Rivelatore Layer1(4, 0.02, 27, Theta_Multi, 0.012, 0.003);
    Rivelatore Layer2(7, 0.02, 27, Theta_Multi, 0.012, 0.003);
    

    //Definiamo una struct che contenga le caratteristiche del vertice (coordinate e molteplicita') 
    typedef struct{
        double x, y, z;
        int molt;
    } Vertice;
    
    static Vertice inizio;
    
    tree->Branch("VertMult", &inizio.x, "x/D:y:z:molt/I"); 
    tree->Branch("Hit1", &riv_1);
    tree->Branch("Hit2", &riv_2);
  	

    tree->SetAutoSave(0); //Rimuove backup cycle
    

    //Creiamo una particella (fuori dal for sul numero degli esperimenti, cosi viene creata una sola volta)
    Particella* part = new Particella();

    //Creiamo un hit temporaneo
    Punto* hit = new Punto();
    
    //for sul numero di esperimenti
    for(int k=0; k<N_esp; k++){
        //Iniziamo a generare il vertice, ci servono 3 coordinate e la molteplicita
        //Unita di misura della lunghezza = cm
        inizio.x = ptr->Gaus(0.,0.01);
        inizio.y = ptr->Gaus(0.,0.01);
		
		switch(gen_z) {
			case 1:
				inizio.z = ptr->Gaus(0.,5.3);
				break;
			case 2:
				inizio.z = ptr->Rndm()*27. - 13.5;
				break;
			case 3:
				inizio.z = z_fissa;
				break;
			default:
				std::cout << "Scelta non valida. Impostato il settaggio di base: distribuzione gaussiana" << std::endl;
				gen_z = 1;
				inizio.z = ptr->Gaus(0.,5.3);
		}
        
        inizio.molt = (ptr->*rndm_molt)(N);
        /*
        if(k==verboseEvent) {
            ofs << "Vertice in (" << inizio.x << ", "<<inizio.y << ", "<<inizio.z << ") e molteplicita " << inizio.molt << std::endl << std::endl;
        }
	*/
	if (k/1000 == 0){ofs << "(" << inizio.x << ", "<<inizio.y << ", "<<inizio.z << ") e molteplicita " << inizio.molt << std::endl;};
	
        int pos1 = 0;
        int pos2 = 0;

	    //for sul numero di particelle
        for(int i=0; i<inizio.molt; i++) {
            //Generiamo i prodotti nel vertice
            part->SetTheta(ptr->RndTheta());
            part->SetPhi(ptr->Rndm()*2.*pi_greco);
            //if (k == verboseEvent){ofs << "Particella Numero " << i+1 << " : ( " << part->GetTheta() << " , " <<  part->GetPhi() << " )" << std::endl;};
            if (k/1000 == 0){ofs << "Particella Numero " << i+1 << " : ( " << part->GetTheta() << " , " <<  part->GetPhi() << " )" << std::endl;}; 
            
            //Trasporto e multiscattering particella per particella
		
            //BEAM PIPE
            *hit = Beam_Pipe.Hit(Punto(inizio.x, inizio.y, inizio.z), part);
            //if (k == verboseEvent){ofs << "Hit sulla BP della particella " << i << " in (" << hit->GetX() << ", " << hit->GetY() << ", " << hit->GetZ() << ")" << endl;}
            *part = (Beam_Pipe.*rndm_scatt)(part);
            //if (k == verboseEvent){ofs << "Multiscattering della particella " << i << "sulla BP in (" << part->GetTheta() << ", " << part->GetPhi() << ")" << endl;}
            
            //LAYER 1
            *hit = Layer1.Hit(*hit, part);
		
	        //Controlliamo che la z del vertice sia all'interno del rivelatore
            if (TMath::Abs(hit -> GetZ())<=((Layer1.GetH())/2.)){

                //if (k == verboseEvent){ofs << "Hit su L1 della particella " << i << " in (" << hit->GetX() << ", " << hit->GetY() << ", " << hit->GetZ() << ")" << endl;}

                //Immagazziniamo lo smearing (coordinale cilindriche)
                new(hit1[pos1]) Segnale((Layer1.*smearing)(hit, i+1));
                /*
		Segnale *tst = (Segnale*)hit1[pos1];
                if (k == verboseEvent){ofs <<"Particella "<<i<<" su L1) Z , phi = "<<tst->GetZ()<<"; "<<tst->GetPhi()<<std::endl;};
		*/
		    
                //Multiscattering
                *part = (Layer1.*rndm_scatt)(part);
                //if (k == verboseEvent){ofs << "Multiscattering della particella " << i << "su L1 in (" << part->GetTheta() << ", " << part->GetPhi() << ")" << endl;}

                //LAYER 2
                *hit = Layer2.Hit(*hit, part);
		    
	            //Controlliamo che la z del vertice sia all'interno del rivelatore
                if(TMath::Abs(hit -> GetZ())<=((Layer2.GetH())/2.)){
                    //if (k == verboseEvent){ofs << "Hit su L2 della particella " << i << " in (" << hit->GetX() << ", " << hit->GetY() << ", " << hit->GetZ() << ")" << endl;}
                    //Immagazziniamo lo smearing
                    new(hit2[pos2]) Segnale((Layer2.*smearing)(hit, i+1));
                    //Segnale *tst2 = (Segnale*)hit2[pos2];
                    //if (k == verboseEvent){ofs <<"Particella "<<i<<" su L2) Z , phi = "<<tst2->GetZ()<<"; "<<tst2->GetPhi()<<std::endl<< std::endl;};

                    pos2++;
                }
                pos1++;
            }

        } //chiusura del for sul numeri di particelle
	    
	    

        //Generazione dei false hit (uniformi in Z e phi) Z in cm
        for(int i=0; i<N_false_hit; i++) {
            new(hit1[pos1+i]) Segnale(-(Layer1.GetH())/2.+(ptr->Rndm())*(Layer1.GetH()),ptr->Rndm()*2*pi_greco, -(i+1));
            new(hit2[pos2+i]) Segnale(-(Layer2.GetH())/2.+(ptr->Rndm())*(Layer2.GetH()),ptr->Rndm()*2*pi_greco, -(i+1));
        }

        // Debug
        //printf("Entries nel TClonesArray: %d\n",riv_1->GetEntries());
        for (int j=0; j<hit1.GetEntries(); j++){
        Segnale *tst = (Segnale*)hit1[j];
        if (k/1000 == 0){ofs <<"Particella "<<j+1<<") Z , phi = "<<tst->GetZ()<<"; "<<tst->GetPhi()<<std::endl;};
        }
        //printf("Entries nel TClonesArray: %d\n",riv_2->GetEntries());
        for (int j=0; j<hit2.GetEntries(); j++){
        Segnale *tst = (Segnale*)hit2[j];
        if (k/1000 == 0){ofs <<"Particella "<<j+1<<") Z , phi = "<<tst->GetZ()<<"; "<<tst->GetPhi()<<std::endl;};
        }
        // fine del debug
        

        tree->Fill();

        riv_1->Clear();
        riv_2->Clear();
    }

    //Salviamo i dati sul file di output 
    Ofile.Write();
    
    //Chiudiamo il file di output
    Ofile.Close();
    
    timer.Stop();
    timer.Print();
	
    //Chiudiamo il file log	
    ofs.close();
	
    //Deallochiamo i puntatori
    delete part;
    delete hit;

}
