#include "TGLViewer.h"
#include "Segnale.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TClonesArray.h"
#include "TStyle.h"
#include "TPolyLine3D.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TGeoManager.h"


void collider_graph(const char* input = "MonteCarlo.root") {
    int dim = 36;
  TClonesArray *riv_1 = new TClonesArray("Segnale",dim);
  TClonesArray *riv_2 = new TClonesArray("Segnale",dim);
  TFile Input_file(input);
  
    // definizione struct
    typedef struct {
        double x, y, z;
        int molt;
    } Vertice;

    static Vertice inizio;

  std::cout << "Evento da studiare:" << std::endl;   
  int N = 0; 
	std::cin >> N;
  
  TObject *obj = (TObject*)Input_file.Get("Generazione");
    
  int N_molt = obj -> GetUniqueID();
  
  TTree *tree = (TTree*)Input_file.Get("Tree");
  TBranch *b1 = tree->GetBranch("VertMult");
  TBranch *b2 = tree->GetBranch("Hit1");
  TBranch *b3 = tree->GetBranch("Hit2");
  b1->SetAddress(&inizio.x);
  b2->SetAddress(&riv_1);
  b3->SetAddress(&riv_2);

  tree->GetEvent(N);
  
  
  
  TCanvas *c1 = new TCanvas();
  
  TGLViewer *view = (TGLViewer*)gPad->GetViewer3D(); //accede al pad del canvas
  TGeoManager *man = new TGeoManager(); //gestisce tutti i volumi presenti
  TGeoHMatrix *tras_rot = new TGeoHMatrix("trans_rot"); //gestisce le rotazioni e le traslazioni ( RotateX(deg) e SetDx(cm) )
  
  tras_rot->RotateY(0.); //non sono sicuro di come pone la terna x,y,z (lo ruoto se penso a z verso l'alto)  87. è bello
  
  //Ora devo creare una BOX che definisca il volume del Master Reference Sistem
  TGeoVolume *top = man->MakeBox("BOX", NULL, 27., 27., 27.); //al posto di Null posso specificare il materiale dell'oggetto, le unità sono in cm
  
  TGeoVolume *beam_pipe = man->MakeTube("PIPE", NULL, 2.96, 3.04, 13.5); //al posto di null potrei creare un medium che registri automaticamente gli hit delle particelle....
  beam_pipe->SetLineColor(kGreen);
  TGeoVolume *riv_1_vol = man->MakeTube("RIV1", NULL, 3.99, 4.01, 13.5);
  riv_1_vol->SetLineColor(kOrange);
  TGeoVolume *riv_2_vol = man->MakeTube("RIV1", NULL, 6.99, 7.01, 13.5);
  riv_2_vol->SetLineColor(kBlue);
  
  man->SetTopVolume(top); //setto il volume madre (esterno)
  top->AddNode(beam_pipe, 0, tras_rot); //creo il layer figlio 0-esimo
  top->AddNode(riv_1_vol, 1, tras_rot);
  top->AddNode(riv_2_vol, 2, tras_rot);
  
  man->CloseGeometry(); //chiudere sempre la geometria creata se no non disegna
  top->Draw("ogl"); //ogl crea l'effettivo oggetto con raytracing ecc, quindi può essere pesante

  vector<TPolyLine3D> tracks; //facilmente utilizzabile per le particelle
    TPolyLine3D pippo = TPolyLine3D();
  for (int j=0;j<inizio.molt;j++){
      tracks.push_back(pippo);
      Segnale* riv_1_data = (Segnale*)riv_1->At(j);
      Segnale* riv_2_data = (Segnale*)riv_2->At(j);
      tracks.at(j).SetLineColor(kRed);
      tracks.at(j).SetPoint(0,inizio.x,inizio.y,inizio.z);
      tracks.at(j).SetPoint(1,4.*TMath::Cos(riv_1_data->GetPhi()),4.*TMath::Sin(riv_1_data->GetPhi()),riv_1_data->GetZ());
      tracks.at(j).SetPoint(2,7.*TMath::Cos(riv_2_data->GetPhi()),7.*TMath::Sin(riv_2_data->GetPhi()),riv_2_data->GetZ());
      tracks.at(j).Draw("same");
  }
}