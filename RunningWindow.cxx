#include "RunningWindow.h"

#include "TMath.h"

using std::vector;

ClassImp(RunningWindow)

bool RunningWindow::raddoppio = 0;

//Copy-Constructor
RunningWindow::RunningWindow(const RunningWindow& source) : TObject(source) {
	//Non allocando memoria possiamo usarne uno semplice
	dmSize=source.dmSize;
	dmStep=source.dmStep;
}

//Operatore =
RunningWindow& RunningWindow::operator=(const RunningWindow& source) {
	//Utilizzo quello che sfrutta il copy constructor per non ripetere tutto
	if(this == &source) return *this;
	this->~RunningWindow();
	new(this) RunningWindow(source);
	return *this;
}


double RunningWindow::media(vector<double> V,int j,double limite) {

    int vec_dim = V.size();
    double temp = 0;
    int count = 0;
    for(int i=j; i<vec_dim && V.at(i)<=limite; i++) {
        temp+=V.at(i);
        count++;
    }
    return temp/(double)count;

}

double RunningWindow::running_window(vector<double> vec,bool &stato_rec) {

    stato_rec = 1; //segna se il vertice e' stato ricostruito o meno
    
    int c_max = 0; //conteggio massimo
    double Z_max = 0; //media delle Zrec nella window con conteggio massimo
    double k_start = 0; //lower bound della finestra
    int j_max = 0; //numero della finestra con conteggio massimo

    if(vec.empty()) {
        stato_rec = 0;
        return 0;
    }
    else {
        double z_0 = vec.at(0)-0.1; //cm //Punto di partenza della prima finestra
        int vec_dim = vec.size();
        for(int j=0; vec.at(vec_dim-1) > z_0 + j*dmStep; j++){

            bool inside = 1; //segna se siamo dentro la finestra
            bool start_window = 1; //ci serve per salvare l'indice del primo elemento del vector che entra nella finestra

            int c = 0; //conteggio
            int k = k_start; //salva l'inizio della finestra

            while((k < vec_dim) && (inside)){
                if((vec.at(k) >= z_0 + j*dmStep) && (vec.at(k) <= z_0 + j*dmStep + dmSize)){
                    if(start_window) {
                        k_start=k; //ci salva da dove partire a scorrere sul vector per la prossima finestra
                        start_window = 0;
                    }
                    c++;
                }
                else if(vec.at(k) > z_0 + j*dmStep + dmSize) inside = 0;
                k++;
            }

            if(c>c_max) {
                c_max = c;
                Z_max = media(vec, k_start, z_0 + j*dmStep + dmSize);
                stato_rec = 1;
                j_max = j;
            }

            //Per ora in caso di parita' viene solo tenuto il caso di finestre uguali vicine, bisogna ancora includere il caso di finestre distanti
            //(basta rimuovere l'if)
            else if(c==c_max) stato_rec = 0;
                /*{
                if(j - j_max == 1 && raddoppio == 0){
                    window = 2 * window;
                    j = 0;
                    j_max = 0;
                    c_max = 0;
                    Z_max = 0;
                    k_start = 0;
                    raddoppio = 1;
                }
                else  stato_rec = 0;
                */
                
            }
            if(stato_rec==0){
                raddoppio = 1;
                dmSize = 2.*dmSize;
                running_window(vec,stato_rec);
            }

        }

        return Z_max;
    }