#include "MyRandom.h"
#include "Riostream.h"
#include "TMath.h"

ClassImp(MyRandom)

bool MyRandom::sFlag = false;

//Costruttore di default
MyRandom::MyRandom() : TRandom3() {
	dmEta = new TH1D();
	dmMolt = new TH1D();
}

//Costruttore standard
MyRandom::MyRandom(const char* input_file, double seed) : TRandom3(seed) {
	TFile F(input_file);
	if(F.IsZombie()) {
		sFlag = true;
	}

	else{
		TH1D* temp_eta = (TH1D*)F.Get("heta");
		TH1D* temp_mol = (TH1D*)F.Get("hmul");
		temp_eta->SetDirectory(0);
		temp_mol->SetDirectory(0);
		F.Close();

		TAxis *xa=temp_eta->GetXaxis();
		Int_t b1=xa->FindBin(-2.);
		Int_t b2=xa->FindBin(2.);
		Double_t xlow=xa->GetBinLowEdge(b1);
		Double_t xhig=xa->GetBinUpEdge(b2);
		Int_t nobins=b2-b1+1;
		dmEta = new TH1D("dmEta","#eta distribution",nobins,xlow,xhig);
		Int_t j=1;
		for(Int_t i=b1;i<=b2;i++)dmEta->SetBinContent(j++,temp_eta->GetBinContent(i));

		dmMolt = temp_mol;
	}
}

//Copy-Constructor
MyRandom::MyRandom(const MyRandom& source) : TRandom3(source) {
	dmEta = new TH1D();
	dmMolt = new TH1D();
	*dmEta = *source.dmEta;
	*dmMolt = *source.dmMolt;
}

//Destructor Non dovrebbe servire deallocare i TH1D* siccome dovrebbe pensarci root
MyRandom::~MyRandom() {}

//Operatore =
MyRandom& MyRandom::operator=(const MyRandom& source) {
    if(this == &source) return *this;
    this->~MyRandom();
    new(this) MyRandom(source);
    return *this;
}

//RndTheta estrae secondo la distribuzione della pseudorapidità e poi calcola theta
double MyRandom::RndTheta() {
	double eta=dmEta->GetRandom();
	double temp=TMath::Exp(-eta);
	return 2.*TMath::ATan(temp);
}
