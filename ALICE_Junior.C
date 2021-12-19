//Librerie custom
#include "MyRandom.h"
#include "Rivelatore.h"
#include "Particella.h"
#include "Vertice.h" //Carica anche Punto.h

void MonteCarlo(usigned int seed = 125) {
    
    MyRandom *ptr = new MyRandom(seed);
    delete gRandom;
    gRandom = ptr;

    //Iniziamo a generare il vertice, ci servono 3 coordinate e la molteplicità
    //Unità di misura della lunghezza = cm
    Vertice Ver(ptr->Gaus(0.,0.01),ptr->Gaus(0.,0.01),ptr->Gaus(0.,5.3),(int) ptr->Rndm());

}