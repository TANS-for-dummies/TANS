//Debugging file, scritto per controllare che lle varie 
//componenti della simulazione facciano il loro lavoro
//Crea un file .root con diversi grafici per eseguire controlli

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

void DebugALICE(bool scat = 1) {
    std::cout << "Inizio dei controlli" << std::endl << std::endl;
    std::cout << "Partiamo verificando che il vertice venga generato bene" << std::endl;
    std::cout << "Controlliamo sia le coordinate che le tre opzioni sulla molteplicità'" << std::endl << std::endl;

    //Costanti
    double pi_greco = TMath::Pi();

    //Settaggi input e output
    unsigned int seed = 0;
    const char* input_file = "kinem.root";
    const char* output_file = "DebugLog.root";
    int NdebugEv = 1000000;

    //Generatore di numeri random
    MyRandom *ptr = new MyRandom(input_file,seed);
    delete gRandom;
    gRandom = ptr;
    if (ptr->GetFlag()) {
    	std::cout << "File con le distribuzioni non trovato" << std::endl; 
    	return;
    }

    //Apriamo il file di output
    TFile Ofile(output_file, "RECREATE");   

    //Rivelatori
    double Theta_Multi = (1/TMath::Sqrt2())*0.001; //rad
	Rivelatore Beam_Pipe(3, 0.08,52, Theta_Multi); //H=52 per contenere tutte le particelle generate con l'accettanza data
	Rivelatore Layer1(4, 0.02, 27, Theta_Multi);
	Rivelatore Layer2(7, 0.02, 27, Theta_Multi); 

    //Funtore per il multiscattering (per controllare che funzioni)
    Particella (Rivelatore::*rndm_scatt) (Particella*, MyRandom*);
    if (scat){rndm_scatt = &Rivelatore::MultiScattering;}
    else {rndm_scatt = &Rivelatore::ZeroScattering;}

    //Iniziamo a controllare la generazione del vertice, guardiamo le distribuzioni
    TH1D *h1 = new TH1D("Vert1","X del Vertice",100,-0.1,0.1);
    TH1D *h2 = new TH1D("Vert2","Y del Vertice",100,-0.1,0.1);
    TH1D *h3 = new TH1D("Vert3","Z del Vertice",100,-25,25);
    TH1D *h4 = new TH1D("Vert4","Molteplicita' costante",101,-0.5,100.5);
    TH1D *h5 = new TH1D("Vert5","Molteplicita' uniforme",101,-0.5,100.5);
    TH1D *h6 = new TH1D("Vert6","Molteplicita' secondo distribuzione",101,-0.5,100.5);
    for(int i=0;i<NdebugEv;i++) {
        h1->Fill(ptr->Gaus(0.,0.01));
        h2->Fill(ptr->Gaus(0.,0.01));
        h3->Fill(ptr->Gaus(0.,5.3));
        h4->Fill(ptr->RndMolt_fissa(100));
        h5->Fill(ptr->RndMolt_unif(100));
        h6->Fill(ptr->RndMolt(100));
    }

    std::cout << "Ora vediamo il percorso della particella" << std::endl << std::endl;
    TH1D *h7 = new TH1D("Part1","Theta iniziale della particella",100,0,3.15);
    TH1D *h8 = new TH1D("Part2","Phi iniziale della particella",100,0.,6.3);
    TH1D *h9 = new TH1D("Eta","Z del primo hit",100,-3.,3);
    TH1D *h10 = new TH1D("Hit1X","X del primo hit",100,0.,4.);
    TH1D *h11 = new TH1D("Hit1Y","Y del primo hit",100,0.,4.);
    TH1D *h12 = new TH1D("Hit1Z","Z del primo hit",100,-15.,15.);
    TH1D *h13 = new TH1D("MultScatterBPtheta","Theta della particella dopo il multiscattering su BP",100,0.,3.15);
    TH1D *h14 = new TH1D("MultScatterBPphi","Phi della particella dopo il multiscattering su BP",100,0.,6.3);
    TH1D *h15 = new TH1D("Hit2X","X del secondo hit",125,0.,5.);
    TH1D *h16 = new TH1D("Hit2Y","Y del secondo hit",125,0.,5.);
    TH1D *h17 = new TH1D("Hit2Z","Z del secondo hit",200,-30.,30.);
    TH1D *h18 = new TH1D("Smear1Z","Z del primo rivelatore dopo lo smearing",200,-30.,30.);
    TH1D *h19 = new TH1D("Smear1Phi","Z del primo rivelatore dopo lo smearing",100,0.,6.3);
    TH1D *h20 = new TH1D("MultScatterL1theta","Theta della particella dopo il multiscattering su L1",100,0.,3.15);
    TH1D *h21 = new TH1D("MultScatterL1phi","Phi della particella dopo il multiscattering su L1",100,0.,6.3);
    TH1D *h22 = new TH1D("Hit3X","X del terzo hit",200,0.,8);
    TH1D *h23 = new TH1D("Hit3Y","Y del terzo hit",200,0.,8);
    TH1D *h24 = new TH1D("Hit3Z","Z del terzo hit",300,-45.,45.);
    TH1D *h25 = new TH1D("Smear2Z","Z del secondo rivelatore dopo lo smearing",300,-45.,45.);
    TH1D *h26 = new TH1D("Smear2Phi","Z del secondo rivelatore dopo lo smearing",100,0.,6.3);

    Punto Vertice(0.,0.,0.);
    Particella* temp_part = new Particella();
    Punto* temp_hit = new Punto();
    Segnale* temp_segnale =new Segnale();

    for(int i=0;i<NdebugEv;i++) {
        //Generiamo la particella con direzione casuale
        temp_part->SetTheta(ptr->RndTheta());
        temp_part->SetPhi(ptr->Rndm()*2.*pi_greco);
        double theta = temp_part->GetTheta();
        double temp = TMath::Tan(theta/2.);
        double eta = -1.*TMath::Log(temp);
        h7->Fill(theta);
        h8->Fill(temp_part->GetPhi());
        h9->Fill(eta);

        //Calcoliamo l'hit sulla beam pipe
        *temp_hit = Beam_Pipe.Hit(Vertice, temp_part);
        h10->Fill(temp_hit->GetX());
        h11->Fill(temp_hit->GetY());
        h12->Fill(temp_hit->GetZ());

        //Ora facciamo il multiscattering sulla beam pipe
        *temp_part=(Beam_Pipe.*rndm_scatt)(temp_part, ptr);
        h13->Fill(temp_part->GetTheta());
        h14->Fill(temp_part->GetPhi());

        //Passiamo ora a calcolare l'hit sul primo rivelatore
        *temp_hit = Layer1.Hit(*temp_hit, temp_part);
        h15->Fill(temp_hit->GetX());
        h16->Fill(temp_hit->GetY());
        h17->Fill(temp_hit->GetZ());
        //Per ora non ci preoccupiamo che Z sia entro le dimensioni del rivelatore
        //Registriamo ogni singolo dato, vogliamo solo verificare che le funzioni
        //si comportino bene

        //Ora facciamo smearing sull'hit per salvare l'output del rivelatore
        *temp_segnale = Segnale(Layer2.Smearing(temp_hit, ptr, i+1));
        h18->Fill(temp_segnale->GetZ());
        h19->Fill(temp_segnale->GetPhi());

        //Ora calcoliamo il multiscattering sul primo rivelatore
        *temp_part = (Layer1.*rndm_scatt)(temp_part, ptr);
        h20->Fill(temp_part->GetTheta());
        h21->Fill(temp_part->GetPhi());

        //Calcoliamo l'hit sul secondo rivelatore
        *temp_hit = Layer2.Hit(*temp_hit, temp_part);
        h22->Fill(temp_hit->GetX());
        h23->Fill(temp_hit->GetY());
        h24->Fill(temp_hit->GetZ());

        //Ora facciamo smearing sull'hit per salvare l'output del rivelatore
        *temp_segnale = Segnale(Layer2.Smearing(temp_hit, ptr, i+1));
        h25->Fill(temp_segnale->GetZ());
        h26->Fill(temp_segnale->GetPhi());
    }    

    std::cout << "Fine dei controlli" << std::endl;

    Ofile.Write();
    Ofile.Close();
}