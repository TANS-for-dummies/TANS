//Librerie custom
#include "Segnale.h"
#include "Tracklet.h"

//Librerie
#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "vector"

using std::vector;

void Ricostruzione_Vertice(int dim = 36){
    //dim = 36: dimensione minima dei TClonesArray

    //Costanti
    double pi_greco = TMath::Pi();

    //Settaggi
    double r1 = 4.; //cm
    double r2 = 7.; //cm
    double delta_phi = 0.4; //ampiezza angolare in rad entro cui cercare corrispondenza hit

    // definizione struct
    typedef struct {
        double x, y, z;
        int molt;
    } Vertice;

    static Vertice inizio;

    // Dichiarazione TClonesArray
    TClonesArray *riv_1 = new TClonesArray("Segnale",dim);
    TClonesArray *riv_2 = new TClonesArray("Segnale",dim);
 
    //Apertura file di input
    TFile Input_file("MonteCarlo.root");

    //Lettura TTree e branch
    TTree *tree = (TTree*)Input_file.Get("Tree");
    TBranch *b1 = tree->GetBranch("VertMult");
    TBranch *b2 = tree->GetBranch("Hit1");
    TBranch *b3 = tree->GetBranch("Hit2");

    // Definizione degli indirizzi per la lettura dei dati su ttree
    b1->SetAddress(&inizio.x);
    b2->SetAddress(&riv_1);
    b3->SetAddress(&riv_2);

    //Creiamo un tracklet, con r1 ed r2 fissati
    Tracklet* tr = new Tracklet(r1, r2, 0., 0.); //Primo punto: layer 1; Secondo punto: layer 2

    //Creiamo un istogramma e un vector per le z
    TH1D* histo_z = new TH1D("histo_z", "Istogramma delle z", 250, -21.2, 21.2); //4 sigma a destra e a sinistra
    vector<double> vec_z; 

    // loop sugli ingressi nel TTree
    for(int i=0; i<tree->GetEntries(); i++){
        tree->GetEvent(i);
        /*
        std::cout << "Evento " << i+1 << "; Molteplicita= " << inizio.molt << std::endl;
        std::cout << "X,Y,Z = " << inizio.x << "; " << inizio.y << "; " << inizio.z << std::endl;

        int num = riv_1->GetEntries();
        std::cout << "Numero di elementi nel primo TClonesArray " << num << std::endl;
        for (int j=0; j<num; j++){
            Segnale *tst=(Segnale*)riv_1->At(j);
            std::cout << "Segnale " << j+1 << ") z, phi = " << tst->GetZ() << "; " << tst->GetPhi() << std::endl;
        }

        num = riv_2->GetEntries();
        std::cout << "Numero di elementi nel secondo TClonesArray " << num << std::endl;
        for (int j=0; j<num; j++){
            Segnale *tst=(Segnale*)riv_2->At(j);
            std::cout << "Segnale " << j+1 << ") z, phi = " << tst->GetZ() << "; " << tst->GetPhi() << std::endl;
        }
        */

       for(int j=0; j<riv_1->GetEntries(); j++){ //for sul layer 1
            Segnale* interazione1 = (Segnale*)riv_1->At(j);
            tr->SetZ1(interazione1->GetZ());

            for(int k=0; k<riv_2->GetEntries(); k++){  //for sul layer 2
                Segnale* interazione2 = (Segnale*)riv_2->At(k);

                if( (interazione1->GetPhi() - interazione2->GetPhi()) <= delta_phi ){
                    tr->SetZ2(interazione2->GetZ());
                    histo_z->Fill(tr->Intersezione());
                    vec_z.push_back(tr->Intersezione()); //riempiamo il vector
                }
           }
       }
        
        histo_z->DrawCopy();
        vec_z.sort(vec_z.begin(), vec_z.end()); //riordiniamo gli elementi del vector in ordine crescente
        
        int max_bin = histo_z->GetMaximumBin(); //bin con massimo numero di conteggi 
        double min_edge = histo_z->GetBinLowEdge(max_bin); //estremo inferiore del bin considerato 
        double max_edge = min_edge + histo_z->GetBinWidth(max_bin); //estremo superiore del bin considerato 
        
        //Calcoliamo la media degli elementi presenti nel bin con massimo numero di conteggi
        int count = 0;
        double media = 0;
        for(int j=0; j<vec_z.size(); j++){
            if(vec_z.at(j)>min_edge && vec_z.at(j)<max_edge) {
                media = media + vec_z.at(j);
                count++;
            }
        }
        
        if(count != 0) media = media/count; //z ricostruita
        

        
       //Reset dell'istogramma e clear del vector
       histo_z->Reset();
       vec_z.clear(); 
    }

}
