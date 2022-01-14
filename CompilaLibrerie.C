#include "TString.h"
#include "TSystem.h"

void CompilaLibrerie(TString Opzione="fast"){
	TString Scelta;
	if(Opzione.Contains("force")){
		Scelta = "kfg";
	}
	else{
		Scelta = "kg";
	}
	gSystem->CompileMacro("Punto.cxx",Scelta.Data());
	gSystem->CompileMacro("Particella.cxx",Scelta.Data());
	gSystem->CompileMacro("Rivelatore.cxx",Scelta.Data());
	gSystem->CompileMacro("MyRandom.cxx",Scelta.Data());
	gSystem->CompileMacro("Segnale.cxx",Scelta.Data());
	gSystem->CompileMacro("Tracklet.cxx",Scelta.Data());
	
	gSystem->CompileMacro("ALICE_Junior.C",Scelta.Data());
	gSystem->CompileMacro("Ricostruzione_Vertice.C",Scelta.Data());
}
