//Librerie custom
#include "Segnale.h"
#include "Tracklet.h"

//Librerie
#include "Riostream.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TGraphErrors.h"
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
    
    const int dim_molt = 10; //numero di molteplicita studiate
    double molteplicita_studiate[] = {3,5,7,9,11,15,20,30,40,50}; //tiene conto delle molteplicita che vogliamo analizzare con i grafici

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
    
    double Z_rec;

    //Creiamo gli istogrammi con cui analizzeremo l'efficienza del nostro algoritmo
    TH1D* deltaZ = new TH1D("deltaZ","Residui",100,-1000,1000);
    deltaZ->GetXaxis()->SetTitle("Zrec-Zvera [#mum]");
    deltaZ->SetMarkerStyle(22);
    deltaZ->SetMarkerColor(52);
    deltaZ->SetLineColor(kBlack);
    
    //creiamo un array di istogrammi
    TH1D *histo_molt[dim_molt];

    char nome[30];
    char titolo[50];

    //loop sull'array di molteplicita studiate per creare gli istogrammi di deltaZ per
    //singole molteplicita
    for (int i=0;i<dim_molt;i++) {
        sprintf(nome,"fixed molt %f",molteplicita_studiate[i]);
        sprintf(titolo,"Residui - molteplicita' fissata a %f",molteplicita_studiate[i]);
        histo_molt[i] = new TH1D(nome,titolo,100,-1000,1000);
        histo_molt[i]->GetXaxis()->SetTitle("Zrec-Zvera [#mum]");
        histo_molt[i]->SetMarkerStyle(33);
    }

    //Creiamo il grafico dell'efficienza
    double s_molt[dim_molt] = {0.}; //array di errori per la molteplicita
    double eff[dim_molt] = {0.};
    double conta_molt[dim_molt] = {0.};//conta quanti eventi hanno una certa molteplicita fissata
    double s_eff[dim_molt] = {0.}; //array di errori per l'efficienza
    TGraphErrors *efficienza;



    
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
        
        //histo_z->DrawCopy();
        
        int max_bin = histo_z->GetMaximumBin(); //bin con massimo numero di conteggi 
        double min_edge = histo_z->GetBinLowEdge(max_bin); //estremo inferiore del bin considerato 
        double max_edge = min_edge + histo_z->GetBinWidth(max_bin); //estremo superiore del bin considerato 
        
        //Calcoliamo una media dei conteggi nei bin non nulli, la useremo per fare un controllo sul vertice
        int count_hist = 0;
        double media_hist = 0.;
        for(int j=0; j<histo_z->GetNbinsX(); j++){
            if(histo_z->GetBinContent(j+1)>0) {
                media_hist = media_hist + histo_z->GetBinContent(j+1);
                count_hist++;
            }
        }
        
        if(count_hist != 0) media_hist = media_hist/count_hist;
        
        int vec_dim = vec_z.size();

        bool Rec = 0;

        if(media_hist<histo_z->GetBinContent(max_bin)){

            //Calcoliamo la media degli elementi presenti nel bin con massimo numero di conteggi
            int count = 0;
            double somma = 0;
            for(int j=0; j<vec_dim; j++){
                if(vec_z.at(j)>min_edge && vec_z.at(j)<max_edge) {
                    somma = somma + vec_z.at(j);
                    count++;
                }
            }
        
            if(count != 0) Z_rec = somma/count;  
            Rec = 1;      
            //std::cout << "Z rec: " << Z_rec << std::endl;  

        }
        else if(vec_dim==1){
            Z_rec=vec_z.at(0);
            Rec = 1; 
            //std::cout << "Z rec: " << Z_rec << std::endl;

        }
            
        if(Rec) {
            deltaZ->Fill((Z_rec-inizio.z)*10000);
            
            for(int j=0;j<dim_molt;j++)
            if(inizio.molt==molteplicita_studiate[j]) {
                histo_molt[j]->Fill((Z_rec-inizio.z)*10000);
                eff[j]++;
            }
        
        }

        for(int j=0;j<dim_molt;j++) {
            if(inizio.molt==molteplicita_studiate[j]) conta_molt[j]++;
        }

       //Reset dell'istogramma e clear del vector
       histo_z->Reset();
       vec_z.clear(); 
    }

    for(int i=0;i<dim_molt;i++) {
        eff[i]=eff[i]/conta_molt[i];
    }

    deltaZ->DrawCopy("pe");

    TCanvas* c2 = new TCanvas("c2","c2",80,80,1500,1000);
    c2->Divide(5,2);
    for(int i=0;i<dim_molt;i++) {
        c2->cd(i+1);
        histo_molt[i]->DrawCopy("pe");
    }

    TCanvas* c3 = new TCanvas("c3","c3",80,80,775,500);
    efficienza = new TGraphErrors(dim_molt,molteplicita_studiate,eff,s_molt,s_eff);
    efficienza->SetTitle("Efficienza vs Molteplicita'");
    efficienza->GetXaxis()->SetTitle("Molteplicita'");
    efficienza->GetYaxis()->SetTitle("Efficienza");
    efficienza->SetMarkerStyle(20);
    efficienza->SetMarkerColor(97);
    efficienza->Draw("APC");
}
