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
#include "TStopwatch.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "vector"
#include "algorithm"
#include "TF1.h"

using std::vector;

bool running_window(vector<double>, double, double&);
double media(vector<double>,int,double);

void Ricostruzione_Vertice(const char* input = 'MonteCarlo.root', double window = 0.5, int n_sigma = 3){ 
    //window in cm
    //n_sigma: numero di deviazioni standard considerate per il taglio sulla Z
    //input: nome del file in input (solo .root)
    
    
    //Costanti
    double pi_greco = TMath::Pi();
    

    //Settaggi
    double r1 = 4.; //cm
    double r2 = 7.; //cm
    double delta_phi = 0.004; //ampiezza angolare in rad entro cui cercare corrispondenza hit
    double sigma_Z = 5.3; //cm (deviazione standard della generazione delle Z

    //Avviamo il timer	
    TStopwatch timer;
    timer.Start();

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
    TFile Input_file(input);

    
    //------------------------------------------------------------Lettura "Generazione"------------------------------------------------------------------
 
    TObject *obj = (TObject*)Input_file.Get("Generazione");
    vector<double> molteplicita_studiate; //è dichiarato double per fare il grafico dopo
    double molteplicita_studiate_standard[10]= {3,5,7,9,11,15,20,30,40,50};
    
    //distribuzione data
    if(obj -> GetUniqueID() == 0){        
        for(int i = 0; i < 10; i++){
            molteplicita_studiate.push_back(molteplicita_studiate_standard[i]); //tiene conto delle molteplicita' medie che vogliamo analizzare con i grafici
        }
    }
    
    //molteplicità fissata
    else if(obj -> GetUniqueID() < 0){
        molteplicita_studiate.push_back(obj->GetUniqueID());
    }
    
    //molteplicità uniforme
    else{
        for(int i = 0; i < 10; i++){
            if(molteplicita_studiate_standard[i] < obj->GetUniqueID()){
                molteplicita_studiate.push_back(molteplicita_studiate_standard[i]);
            }
            molteplicita_studiate.push_back(obj->GetUniqueID());
        }
    }
    
    const int dim_molt = molteplicita_studiate.size();
    
    //----------------------------------------------------------------------------------------------------------------------------------------------------
    
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
    

    //Creiamo gli istogrammi con cui analizzeremo la risoluzione del nostro algoritmo
    TH1D* deltaZ = new TH1D("deltaZ","Residui",200,-1000,1000);
    deltaZ->GetXaxis()->SetTitle("Zrec-Zvera [µm]");
    deltaZ->SetMarkerStyle(22);
    deltaZ->SetMarkerColor(52);
    deltaZ->SetLineColor(kBlack);
    
    //creiamo un array di istogrammi
    TH1D *histo_molt[dim_molt];

    char nome[30];
    char titolo[80];

    //loop sull'array di molteplicita' studiate per creare gli istogrammi di deltaZ per
    //singole molteplicita' ed inizializzare le gaussiane
    for (int i=0;i<dim_molt; i++) {
        sprintf(nome, "fixed molt center %f", molteplicita_studiate.at(i));
        sprintf(titolo,"Residui - molteplicita' fissata da %f a %f", molteplicita_studiate.at(i) - 0.5, molteplicita_studiate.at(i) - 0.5);
        histo_molt[i] = new TH1D(nome, titolo, 400, -1000, 1000);
        histo_molt[i] -> GetXaxis() -> SetTitle("Zrec-Zvera [#mum]");
        histo_molt[i] -> SetMarkerStyle(33);
    }
    

    //Creiamo il grafico dell'efficienza e della risoluzione in funzione della molteplicita'
    double s_molt[dim_molt] = {0.}; //array di errori per la molteplicita
    double eff[dim_molt] = {0.};
    double conta_molt[dim_molt] = {0.};//conta quanti eventi hanno una certa molteplicita fissata
    double s_eff[dim_molt] = {0.}; //array di errori per l'efficienza
    double ris[dim_molt] = {0.}; //array per le risoluzioni prese dai fit
    double s_ris[dim_molt] = {0.}; //array di errori sulla risoluzione 
    TGraphErrors *efficienza;
    TGraphErrors *risoluzione;


    //Creiamo il grafico dell'efficienza e della risoluzione in funzione di Z
    const int dim_Z = 15;
    double Z[dim_Z] = {-12.6, -10.8, -9., -7.2, -5.4, -3.6, -1.8, 0., 1.8, 3.6, 5.4, 7.2, 9., 10.8, 12.6};
    double s_Z[dim_Z] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}; 
    double eff_Z[dim_Z] = {0.};
    double conta_Z[dim_Z] = {0.};//conta quanti eventi ci sono in un certo intervallo di Z
    double s_eff_Z[dim_Z] = {0.}; //array di errori per l'efficienza
    double ris_Z[dim_Z] = {0.}; //array per le risoluzioni prese dai fit
    double s_ris_Z[dim_Z] = {0.}; //array di errori sulla risoluzione 
    TGraphErrors *efficienza_Z;
    TGraphErrors *risoluzione_Z;



    
    // loop sugli ingressi nel TTree
    for(int i=0; i<tree->GetEntries(); i++){
        tree->GetEvent(i);
        
        //Controlliamo che le Z generate siano entro n_sigma
        if(TMath::Abs(inizio.z) < (n_sigma*sigma_Z)){
        
           for(int j=0; j<riv_1->GetEntries(); j++){ //for sul layer 1
                Segnale* interazione1 = (Segnale*)riv_1->At(j);
                tr->SetZ1(interazione1->GetZ());

                for(int k=0; k<riv_2->GetEntries(); k++){  //for sul layer 2
                    Segnale* interazione2 = (Segnale*)riv_2->At(k);

                    if( TMath::Abs(interazione1->GetPhi() - interazione2->GetPhi()) <= delta_phi ){
                        tr->SetZ2(interazione2->GetZ());
                        histo_z->Fill(tr->Intersezione());//Riempiamo l'istogramma con le intersezioni tra tracklet ed asse del fascio
                        vec_z.push_back(tr->Intersezione()); //riempiamo il vector
                    }
               }
           } //fine del loop sul TClonesArray

            //histo_z->DrawCopy();

            sort(vec_z.begin(), vec_z.end()); //riordiniamo il vector in ordine crescente

            bool Rec = 1; //indica che riusciamo a ricostruire il vertice
            double Z_rec = 0;

            Rec = running_window(vec_z, window, Z_rec); //Ricostruzione con metodo della running window

            if(Rec) {deltaZ -> Fill((Z_rec-inizio.z)*10000);}


            //DA COMPLETARE (?)
            for(int j=0; j<TMath::Max(dim_molt, dim_Z); j++){
                if(j<dim_molt && (inizio.molt>molt_min[j]) && (inizio.molt<molt_max[j])) {
                    conta_molt[j]++;
                    if(Rec){
                        histo_molt[j]->Fill((Z_rec-inizio.z)*10000);
                        eff[j]++;
                    }
                }
                
                if(j<dim_Z && (inizio.z>(Z[j]-s_Z[j])) && (inizio.z<(Z[j]+s_Z[j]))) {
                    conta_Z[j]++;
                    if(Rec){
                        //histo_molt[j]->Fill((Z_rec-inizio.z)*10000);
                        eff_Z[j]++;
                    }
                }
            }
        
           //Reset dell'istogramma e clear del vector
           histo_z->Reset();
           vec_z.clear(); 
        }  //chiusura if sulle z
    } //chiusura del for sugli eventi

    deltaZ->DrawCopy("pe");

    
    TCanvas* c2 = new TCanvas("c2","c2",80,80,1500,1000);
    c2->Divide(5,2);

    for(int i=0; i<TMath::Max(dim_molt, dim_Z); i++) {
        
        if(i<dim_molt){
        
            //Calcolo efficienza e relativo errore
            eff[i]=eff[i]/conta_molt[i];
            s_molt[i] = 0.5;
            double s_eff[i] = TMath::Sqrt((1.-eff[i])*eff[i]/conta_molt[i]);//Errore binomiale


            //Fit delle gaussiane
            c2->cd(i+1);
            histo_molt[i] -> Fit("gaus");
            TF1 *Gauss = histo_molt[i] -> GetFunction("gaus");
            gStyle->SetOptFit(1111);
            histo_molt[i]->DrawCopy("pe");


            //Calcolo della risoluzione
            ris[i] = Gauss -> GetParameter(2);
            s_ris[i] = Gauss -> GetParError(2);
            
        }    
        
        if(i<dim_Z) {
            
            //Calcolo efficienza e relativo errore
            eff_Z[i]=eff_Z[i]/conta_Z[i];
            double s_eff_Z[i] = TMath::Sqrt((1.-eff_Z[i])*eff_Z[i]/conta_Z[i]);//Errore binomiale
        }   

    }

    
    //Efficienza in funzione della molteplicità
    TCanvas* c3 = new TCanvas("c3","c3",80,80,775,500);
    efficienza = new TGraphErrors(dim_molt,&(molteplicita_studiate[0]),eff,s_molt,s_eff);
    efficienza->SetTitle("Efficienza vs Molteplicita'");
    efficienza->GetXaxis()->SetTitle("Molteplicita'");
    efficienza->GetYaxis()->SetTitle("Efficienza");
    efficienza->SetMarkerStyle(20);
    efficienza->SetMarkerColor(97);
    efficienza->Draw("APC");
    
    
    //Efficienza in funzione di Z
    TCanvas* c4 = new TCanvas("c4","c4",80,80,775,500);
    efficienza_Z = new TGraphErrors(dim_Z, Z, eff_Z, s_Z, s_eff_Z);
    efficienza_Z->SetTitle("Efficienza vs Z");
    efficienza_Z->GetXaxis()->SetTitle("Z");
    efficienza_Z->GetYaxis()->SetTitle("Efficienza");
    efficienza_Z->SetMarkerStyle(20);
    efficienza_Z->SetMarkerColor(97);
    efficienza_Z->Draw("APC");
    
    
    //Risoluzione
    TCanvas* c5 = new TCanvas("c5","c5",80,80,775,500);
    risoluzione = new TGraphErrors(dim_molt,&(molteplicita_studiate[0]),ris,s_molt,s_ris);
    risoluzione->SetTitle("Risoluzione vs Molteplicita'");
    risoluzione->GetXaxis()->SetTitle("Molteplicita'");
    risoluzione->GetYaxis()->SetTitle("Risoluzione (µm)");
    risoluzione->SetMarkerStyle(33);
    risoluzione->SetMarkerColor(77);
    risoluzione ->Draw ("APC");



    timer.Stop();
    timer.Print();
}












bool running_window(vector<double> vec,double window,double &Z) {

    bool stato_rec = 1; //segna se il vertice è stato ricostruito o meno
    double step = window/2.; //di quanto si sposta la finestra ad ogni incremento
    int c_max = 0; //conteggio massimo
    double Z_max = 0; //media delle Zrec nella window con conteggio massimo
    double k_start = 0; //lower bound della finestra
    int j_max = 0; //numero della finestra con conteggio massimo
    bool raddoppio = 0; //segna se è già stata raddoppiata la finestra una volta

    if(vec.empty()) return 0;
    else {
        double z_0 = vec.at(0);
        int vec_dim = vec.size();
        for(int j=0; vec.at(vec_dim-1) > z_0 + j*step; j++){

            bool inside = 1; //segna se siamo dentro la finestra
            bool start_window = 1; //ci serve per salvare l'indice del primo elemento del vector che entra nella finestra
            int c = 0; //conteggio
            int k = k_start; //salvaq l'inizio della finestra

            while((k < vec_dim) && (inside)){
                if((vec.at(k) >= z_0 + j*step) && (vec.at(k) <= z_0 + j*step + window)){
                    if(start_window) {
                        k_start=k; //ci salva da dove partire a scorrere sul vector per la prossima finestra
                        start_window = 0;
                    }
                    c++;
                }
                else if(vec.at(k) > z_0 + j*step + window) inside = 0;
                k++;
            }

            if(c>c_max) {
                c_max = c;
                Z_max = media(vec, k_start, z_0 + j*step + window);
                stato_rec = 1;
                j_max = j;
            }

            //Per ora in caso di parità viene solo tenuto il caso di finestre uguali vicine, bisogna ancora includere il caso di finestre distanti
            //(basta rimuovere l'if)
            else if(c==c_max){
                if(j - j_max = 1 && raddoppio = 0){
                    window = 2 * window;
                    j = 0;
                    j_max = 0;
                    c_max = 0;
                    z_max = 0;
                    k_start = 0;
                    raddoppio = 1;
                 else  stato_rec = 0;
                }
                
            }

        }

        Z=Z_max;
        return stato_rec;
    }
}







double media(vector<double> V,int j,double limite) {

    int vec_dim = V.size();
    double temp = 0;
    int count = 0;
    for(int i=j; i<vec_dim && V.at(i)<=limite; i++) {
        temp+=V.at(i);
        count++;
    }
    return temp/(double)count;

}
