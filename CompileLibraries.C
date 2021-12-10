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
	gSystem->CompileMacro("Vertice.cxx",Scelta.Data());
}