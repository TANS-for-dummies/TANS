#include "Rivelatore.h"
#include "TMath.h"


ClassImp(Rivelatore)

//Costruttore di default
Rivelatore::Rivelatore(): dmR(0.), dmS(0.), dmH(0.), dmTheta(1), TObject() {}

//Costruttore standard
Rivelatore::Rivelatore(double r, double s, double H, double Theta): dmR(r), dmS(s), dmH(H), dmTheta(Theta), TObject() {}

//Copy
Rivelatore::Rivelatore(const Rivelatore& source) : dmR(source.dmR), dmS(source.dmS), dmH(source.dmH), dmTheta(source.dmTheta), TObject(source) {}

//Distruttore
Rivelatore::~Rivelatore() {}

//Operatore copia
Rivelatore& Rivelatore::operator=(const Rivelatore& source){
    if(this == &source) return *this;
    this->~Rivelatore();
    new(this) Rivelatore(source);
    return *this;
}

Particella Rivelatore::MultiScattering(Particella *part, MyRandom *ptr){
    //La particella è considerata descritta da 2 angoli
    double ThetaP = ptr -> Gaus(dmTheta,(TMath::Sqrt2())*dmTheta); //DA CHIEDERE A MASERA CON URGENZA!!!!!!!!!!!!!!!!!
    double PhiP = ptr -> Rndm()*2.*TMath::Pi();

    double mr[3][3];
    mr[0][0] = - TMath::Sin(part->GetCoord2()); 
    mr[1][0] = TMath::Cos(part->GetCoord2()); 
    mr[2][0] = 0;
    mr[0][1] = - TMath::Cos(part->GetCoord2())*TMath::Cos(part->GetCoord1()); 
    mr[1][1] = - TMath::Cos(part->GetCoord1())*TMath::Sin(part->GetCoord2());
    mr[2][1] = TMath::Sin(part->GetCoord1()); 
    mr[0][2] = TMath::Cos(part->GetCoord2())*TMath::Sin(part->GetCoord1());
    mr[1][2] = TMath::Sin(part->GetCoord1())*TMath::Sin(part->GetCoord2());
    mr[2][2] = TMath::Cos(part->GetCoord1());

    double scat[3];
    scat[0] = TMath::Sin(ThetaP)*TMath::Cos(PhiP);
    scat[1] = TMath::Sin(ThetaP)*TMath::Sin(PhiP);
    scat[2] = TMath::Cos(ThetaP);

    double final_dir[3];

    for (int i = 0; i < 3; i++){
        final_dir[i] = 0.;
        for (int j = 0; j < 3; j++){
            final_dir[i]+=mr[i][j]*scat[j];
        }
    }
    double final_theta = TMath::ACos(final_dir[2]);
    double final_phi = TMath::ACos(final_dir[0]/(TMath::Sin(final_theta)));
    return Particella(final_theta,final_phi);
}

Particella Rivelatore:: Smearing(Punto *P, MyRandom *ptr){
    //Qui la particella è descritta in coordinate cilindriche, per nostra convenzione (Z,phi)

    double z = P->GetZ() + ptr->Gaus(0,0.012);
    double phi = TMath::ATan( (P->GetY())/(P->GetX()) ) + (ptr->Gaus(0,0.003))/P->GetRadiusXY();

    Particella temp(z,phi);

    return temp; 
}


Punto Rivelatore::Hit(Punto P, Particella *part){
    double x0 = P.GetX();
    double y0 = P.GetY();
    double z0 = P.GetZ();

    double c1 = TMath::Sin(part->GetCoord1()) * TMath::Cos(part->GetCoord2());
    double c2 = TMath::Sin(part->GetCoord1()) * TMath::Sin(part->GetCoord2());
    double c3 = TMath::Cos(part->GetCoord1());

    double r = dmR + 0.5*dmS;

    double delta = (x0*c1 + y0*c2)*(x0*c1 + y0*c2) - (c1*c1 + c2*c2)*(x0*x0 + y0*y0 - r*r);
    double t1 = (-(x0*c1 + y0*c2) + TMath::Sqrt(delta))/(c1*c1 + c2*c2);
    double t2 = (-(x0*c1 + y0*c2) - TMath::Sqrt(delta))/(c1*c1 + c2*c2);
    double t;

    if(t1>0) t=t1;
    else t=t2;

    double x = x0 + c1*t;
    double y = y0 + c2*t;
    double z = z0 + c3*t;

    Punto hit(x,y,z);

    return hit;


}