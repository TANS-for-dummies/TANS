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

bool running_window_1(vector<double>, double, double&);
bool running_window_2(vector<double>, double, double&);
bool rec_hist(TH1D* , vector<double>, double&);
double media(vector<double>,int,double);

void Ricostruzione_Vertice(int dim = 36, double window = 0.5, int n_sigma = 3){ 
    //dim = 36: dimensione minima dei TClonesArray
    //window in cm
    //n_sigma: numero di deviazioni standard considerate per il taglio sulla Z

    //Costanti
    double pi_greco = TMath::Pi();

    //Settaggi
    double r1 = 4.; //cm
    double r2 = 7.; //cm
    double delta_phi = 0.004; //ampiezza angolare in rad entro cui cercare corrispondenza hit
    double sigma_Z = 5.3; //cm (deviazione standard della generazione delle Z)
    
    const int dim_molt = 10; //numero di molteplicita studiate
    double molteplicita_studiate[] = {3,5,7,9,11,15,20,30,40,50}; //tiene conto delle molteplicita' medie che vogliamo analizzare con i grafici
    //double molt_min[] = {2.5,4.5,6.5,8.5,10.5,13.5,17.5,25.5,35.5,45.5};//Minimo e massimo degli intevalli di molteplicita' studiati
    //double molt_max[] = {3.5,5.5,7.5,9.5,11.5,16.5,22.5,34.5,44.5,54.5};
    double molt_min[] = {2.5,4.5,6.5,8.5,10.5,14.5,19.5,29.5,39.5,49.5};//Minimo e massimo degli intevalli di molteplicita' studiati
    double molt_max[] = {3.5,5.5,7.5,9.5,11.5,15.5,20.5,30.5,40.5,50.5};

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
    for (int i=0;i<dim_molt;i++) {
        sprintf(nome,"fixed molt center %f",molteplicita_studiate[i]);
        sprintf(titolo,"Residui - molteplicita' fissata da %f a %f",molt_min[i],molt_max[i]);
        histo_molt[i] = new TH1D(nome,titolo,400,-1000,1000);
        histo_molt[i]->GetXaxis()->SetTitle("Zrec-Zvera [#mum]");
        histo_molt[i]->SetMarkerStyle(33);
    }

    //Creiamo il grafico dell'efficienza e della risoluzione
    double s_molt[dim_molt] = {0.}; //array di errori per la molteplicita
    double eff[dim_molt] = {0.};
    double conta_molt[dim_molt] = {0.};//conta quanti eventi hanno una certa molteplicita fissata
    double s_eff[dim_molt] = {0.}; //array di errori per l'efficienza
    double ris[dim_molt] = {0.}; //array per le risoluzioni prese dai fit
    double s_ris[dim_molt] = {0.}; //array di errori sulla risoluzione 
    TGraphErrors *efficienza;
    TGraphErrors *risoluzione;



    
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

            //Rec=rec_hist(histo_z, vec_z, Z_rec);//Ricostruzione con metodo dell'istogramma
            //Rec=running_window_1(vec_z, window, Z_rec);//Ricostruzione con metodo della running window versione 1
            Rec=running_window_2(vec_z, window, Z_rec);//Ricostruzione con metodo della running window versione 2

            if(Rec) {deltaZ->Fill((Z_rec-inizio.z)*10000);}

            for(int j=0;j<dim_molt;j++){
                if((inizio.molt>molt_min[j])&&(inizio.molt<molt_max[j])) {
                    conta_molt[j]++;
                    if(Rec){
                        histo_molt[j]->Fill((Z_rec-inizio.z)*10000);
                        eff[j]++;
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

    for(int i=0;i<dim_molt;i++) {
        
        //Calcolo efficienza e relativo errore
        eff[i]=eff[i]/conta_molt[i];
        s_molt[i]=(molt_max[i]-molt_min[i])/2.;
        double temp_1 = TMath::Sqrt((1.-eff[i])*eff[i]/conta_molt[i]);//Errore binomiale
        double temp_2 = 1./conta_molt[i];//Errore minimo
        if(temp_1<temp_2)s_eff[i]=temp_2;
        else s_eff[i]=temp_1;


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

    //Efficienza
    TCanvas* c3 = new TCanvas("c3","c3",80,80,775,500);
    efficienza = new TGraphErrors(dim_molt,molteplicita_studiate,eff,s_molt,s_eff);
    efficienza->SetTitle("Efficienza vs Molteplicita'");
    efficienza->GetXaxis()->SetTitle("Molteplicita'");
    efficienza->GetYaxis()->SetTitle("Efficienza");
    efficienza->SetMarkerStyle(20);
    efficienza->SetMarkerColor(97);
    efficienza->Draw("APC");

    //Risoluzione
    TCanvas* c4 = new TCanvas("c4","c4",80,80,775,500);
    risoluzione = new TGraphErrors(dim_molt,molteplicita_studiate,ris,s_molt,s_ris);
    risoluzione->SetTitle("Risoluzione vs Molteplicita'");
    risoluzione->GetXaxis()->SetTitle("Molteplicita'");
    risoluzione->GetYaxis()->SetTitle("Risoluzione (µm)");
    risoluzione->SetMarkerStyle(33);
    risoluzione->SetMarkerColor(77);
    risoluzione ->Draw ("APC");



    timer.Stop();
    timer.Print();
}












bool rec_hist(TH1D *h,vector<double> vec,double Z) {
    
    bool stato_rec=1;

    int max_bin = h->GetMaximumBin(); //bin con massimo numero di conteggi 
    double min_edge = h->GetBinLowEdge(max_bin); //estremo inferiore del bin considerato 
    double max_edge = min_edge + h->GetBinWidth(max_bin); //estremo superiore del bin considerato

    //Calcoliamo una media dei conteggi nei bin non nulli, la useremo per fare un controllo sul vertice
    int count_hist = 0;
    double media_hist = 0.;
    for(int j=0; j<h->GetNbinsX(); j++){
        if(h->GetBinContent(j+1)>0) {
            media_hist = media_hist + h->GetBinContent(j+1);
            count_hist++;
        }
    }

    if(count_hist != 0) media_hist = media_hist/count_hist;
    int vec_dim = vec.size();

    if(media_hist<h->GetBinContent(max_bin)){
        //Calcoliamo la media degli elementi presenti nel bin con massimo numero di conteggi
        int count = 0;
        double somma = 0;
        for(int j=0; j<vec_dim; j++){
            if(vec.at(j)>min_edge && vec.at(j)<max_edge) {
                somma = somma + vec.at(j);
                count++;
            }
        }
        if(count != 0) Z = somma/count;       
    }
    else if(vec_dim==1){
        Z=vec.at(0);
    } 
    else stato_rec = 0;
    
    return stato_rec;
}









bool running_window_1(vector<double> vec,double window,double &Z) {

    bool stato_rec=1;

    //running window: scorriamo il vector cercando la finestra con il maggior numero di conteggi
    int vec_dim = vec.size();
    int j_max = 0; //inizio della finestra contenente il massimo numero di conteggi
    int conteggi_max = 0; //tiene conto dei conteggi nella fienstra contenente il picco
    for(int j=0; j<vec_dim; j++){
        bool inside = 1; //verifica se si e' all'interno della window
        int k = j;
        int conteggi_temp = 0; //tiene conto dei conteggi nella fienstra considerata di volta in volta
        while(inside){
            if( k<vec_dim && (vec.at(k)<=vec.at(j) + window) ){ //controllo di essere dentro la finestra e di non eccedere la dimensione del vector
                conteggi_temp++;
            }
            else inside = 0;
            k++;
        }
        //cout << conteggi_temp << " in " << j << " che sta tra " << vec_z.at(j) << " e " << vec_z.at(j)+window << endl;
        //se troviamo una finestra con piu conteggi del massimo precedente: aggiorniamo conteggi_max e j_max
        if(conteggi_temp > conteggi_max){
            conteggi_max = conteggi_temp;
            j_max = j;
            stato_rec = 1; //tiene conto del fatto che Rec potrebbe essere messo a 0 se in precedenza abbiamo trovato due picchi troppo distanti tra loro
        }
        
        //se troviamo lo stesso numero di conteggi del massimo trovato in precedenza:
        else if(conteggi_temp == conteggi_max){
            //se le due finestre sono distanti:
            if( (vec.at(j) - vec.at(j_max)) >= window ) {
                stato_rec = 0; //non ricostruiamo il vertice
                }
        }
        
    }
    if(stato_rec && vec_dim != 0){
        for(int j=j_max; j<(j_max+conteggi_max); j++){
            Z += vec.at(j)/(double)conteggi_max;
        }
    }   
    else stato_rec = 0;    
    
    return stato_rec;
}








bool running_window_2(vector<double> vec,double window,double &Z) {

    bool stato_rec = 1;
    double step = 0.15;
    int c_max = 0;
    double Z_max = 0; //media delle Zrec nella window con conteggio massimo
    double k_start = 0;

    if(vec.empty()) return 0;
    else {
        double z_0 = vec.at(0);
        int vec_dim = vec.size();
        for(int j=0; vec.at(vec_dim-1) > z_0 + j*step; j++){

            bool inside = 1;
            bool start_window = 1; //ci serve per salvare l'indice del primo elemento del vector che entra nella finestra
            int c = 0;
            int k = k_start;

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
                Z_max = media(vec,k_start,z_0 + j*step + window);
                stato_rec = 1;
            }

            else if(c==c_max){
                double temp_Z = media(vec,k_start,z_0 + j*step + window);
                if(TMath::Abs(Z_max-temp_Z) > 0.15)  stato_rec = 0;
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
