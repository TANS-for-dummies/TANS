#include "Vertice.h"

ClassImp(Vertice)

//Costruttore di default
Vertice::Vertice() : dmX(0.), dmY(0.), dmZ(0.), dmN(0.){}

//Costruttore standard
Vertice::Vertice(double x, double y, double z, int N): dmX(x), dmY(Y), dmZ(z), dmN(N), TObject(){
    dmDir.reserve(N);
}

//Costruttore d copia
Vertice::Vertice(const Vertice& source) : dmX(source.dmX), dmY(source.dmY), dmZ(source.dmZ), dmN(source.dmN), TObject(source){
    if(dmN>0){
        std::vector<Particella> dmDir(dmN);

        for(int i = 0; i<dmN; i++){
            dmDir.At(i)=source.dmDir.At(i);
        }
    }
    else{
        std::vector<Particella> dmDir(dmN);
    }
}

//Distruttore standard
Vertice::~Vertice(){
    delete dmDir;
}

//Operatore copia
Vertice& Vertice::operator=(const Vertice& source){
    if(this == &source) return *this;
    this->~Vertice();
    new(this) Vertice(source);
    return *this;
}

//Setter di dmDir
void Vertice::AddPart(double Theta, double Phi, double p){
    dmDir.push_back(Particella(Theta, Phi, p));
}