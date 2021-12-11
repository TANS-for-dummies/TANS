#include "Vertice.h"

ClassImp(Vertice)

//Costruttore di default
Vertice::Vertice() : Punto(), dmN(0.){}

//Costruttore standard
Vertice::Vertice(double x, double y, double z, int N): Punto(x,y,z), dmN(N){
    dmDir.reserve(N);
}

//Costruttore d copia
Vertice::Vertice(const Vertice& source) : Punto(source), dmN(source.dmN){
    if(dmN>0){
        std::vector<Particella> dmDir(dmN);

        for(int i = 0; i<dmN; i++){
            dmDir.at(i)=source.dmDir.at(i);
        }
    }
    else{
        std::vector<Particella> dmDir(dmN);
    }
}

//Distruttore standard
Vertice::~Vertice(){
    //Per i vector non serve fare il delete
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
