import "TGLViewer.h"



void collider_graph()
{
  TCanvas *c1 = new TCanvas();
  
  TGLViewer *view = (TGLViewer*)gPad->GetViewer3D(); //accede al pad del canvas
  TGeoManager *man = new TGeoManager(); //gestisce tutti i volumi presenti
  TGeoHMatrix *tras_rot = new TGeoHMatrix("trans_rot"); //gestisce le rotazioni e le traslazioni ( RotateX(deg) e SetDx(cm) )
  
  tras_rot->RotateX(90.); //non sono sicuro di come pone la terna x,y,z (lo ruoto se penso a z verso l'alto)
  
  //Ora devo creare una BOX che definisca il volume del Master Reference Sistem
  TGeoVolume *top = man->MakeBox("BOX", NULL, 27., 27., 27.); //al posto di Null posso specificare il materiale dell'oggetto, le unità sono in cm
  
  TGeoVolume *beam_pipe = man->MakeTube("PIPE", NULL, 3., 3.08, 27.); //al posto di null potrei creare un medium che registri automaticamente gli hit delle particelle....
  beam_pipe->SetLineColor(kGreen);
  
  man->SetTopVolume(top); //setto il volume madre (esterno)
  top->AddNode(beam_pipe, 0, trans_rot); //creo il layer figlio 0-esimo
  
  man->CloseGeometry(); //chiudere sempre la geometria creata se no non disegna
  top->Draw("ogl"); //ogl crea l'effettivo oggetto con raytracing ecc, quindi può essere pesante

  TPolyLine3D *l1 = new TPolyLine3D(); //facilmente utilizzabile per le particelle
  l1->SetLineColor(kRed);
  l1->SetPoint(0, 13.5, 0.,0.); //indice, x, y, z
  l1->SetPoint(1, 19.7, 2.78, 0.55);
  l1->SetPoint(2, 25.8, 9.5, 9.5);
  l1->Draw("same");
}
