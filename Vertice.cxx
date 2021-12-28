#include "Vertice.h"

ClassImp(Vertice)

//Costruttore di default
Vertice::Vertice() : Punto(), dmN(0){}

//Costruttore standard
Vertice::Vertice(double x, double y, double z, int N): Punto(x,y,z), dmN(N){
    dmDir.reserve(N); //alloca memoria per il vettore di particelle (grandezza minima allocata nella heap)
}

//Costruttore di copia
Vertice::Vertice(const Vertice& source) : Punto(source), dmN(source.dmN){
    if(dmN>0){
        std::vector<Particella> dmDir(dmN);

        for(int i = 0; i<dmN; i++){
            dmDir.at(i)=source.dmDir.at(i); //at(i) dà la reference dell'i-esimo elemento del vector 
        }
    }
    else{
        std::vector<Particella> dmDir(dmN);
    }
}

//Distruttore standard
Vertice::~Vertice(){
    //Per i vector non serve fare il delete, se ne occupa vector una volta fuori dallo scope
}

//Operatore copia
Vertice& Vertice::operator=(const Vertice& source){
    if(this == &source) return *this;
    this->~Vertice();
    new(this) Vertice(source);
    return *this;
}

//Setter di dmDir
void Vertice::AddPart(double Theta, double Phi){
    dmDir.push_back(Particella(Theta, Phi)); //push_back aggiuge una particella al fondo del vettore dmDir
}
