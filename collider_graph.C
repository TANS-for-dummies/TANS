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


void collider_graph()
{
  TClonesArray *riv_1 = new TClonesArray("Segnale",dim);
  TClonesArray *riv_2 = new TClonesArray("Segnale",dim);
  TFile Input_file(input);
  
  std::cout << "Evento da studiare:" << std::endl;    
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

  tree->GetEvent(i)
  
  
  
  TCanvas *c1 = new TCanvas();
  
  TGLViewer *view = (TGLViewer*)gPad->GetViewer3D(); //accede al pad del canvas
  TGeoManager *man = new TGeoManager(); //gestisce tutti i volumi presenti
  TGeoHMatrix *tras_rot = new TGeoHMatrix("trans_rot"); //gestisce le rotazioni e le traslazioni ( RotateX(deg) e SetDx(cm) )
  
  tras_rot->RotateX(90.); //non sono sicuro di come pone la terna x,y,z (lo ruoto se penso a z verso l'alto)
  
  //Ora devo creare una BOX che definisca il volume del Master Reference Sistem
  TGeoVolume *top = man->MakeBox("BOX", NULL, 27., 27., 27.); //al posto di Null posso specificare il materiale dell'oggetto, le unità sono in cm
  
  TGeoVolume *beam_pipe = man->MakeTube("PIPE", NULL, 3., 3.08, 27.); //al posto di null potrei creare un medium che registri automaticamente gli hit delle particelle....
  beam_pipe->SetLineColor(kGreen);
  TGeoVolume *riv_1 = man->MakeTube("RIV1", NULL, 4., 4.02, 27.);
  riv_1->SetLineColor(kOrange);
  TGeoVolume *riv_2 = man->MakeTube("RIV1", NULL, 7., 7.02, 27.);
  riv_2->SetLineColor(kBlue);
  
  man->SetTopVolume(top); //setto il volume madre (esterno)
  top->AddNode(beam_pipe, 0, trans_rot); //creo il layer figlio 0-esimo
  top->AddNode(riv_1, 1, trans_rot);
  top->AddNode(riv_2, 2, trans_rot);
  
  man->CloseGeometry(); //chiudere sempre la geometria creata se no non disegna
  top->Draw("ogl"); //ogl crea l'effettivo oggetto con raytracing ecc, quindi può essere pesante

  TPolyLine3D *l1 = new TPolyLine3D(); //facilmente utilizzabile per le particelle
  l1->SetLineColor(kRed);
  l1->SetPoint(0, 13.5, 0.,0.); //indice, x, y, z
  l1->SetPoint(1, 19.7, 2.78, 0.55);
  l1->SetPoint(2, 25.8, 9.5, 9.5);
  l1->Draw("same");
}
