#include "MyRandom.h"
#include "Riostream.h"
#include "TMath.h"

ClassImp(MyRandom)

bool MyRandom::sFlag = 0;

//Costruttore di default
MyRandom::MyRandom() : TRandom3() {
	dmEta = new TH1D();
	dmMolt = new TH1D();
}

//Costruttore standard
MyRandom::MyRandom(const char* input_file, double seed) : TRandom3(seed) {
	TFile F(input_file);
	if(F.IsZombie()) {
		sFlag = 1;
	}
	TH1D* temp_eta = (TH1D*)F.Get("heta");
	TH1D* temp_mol = (TH1D*)F.Get("hmul");
	temp_eta->SetDirectory(0);
	temp_mol->SetDirectory(0);
	F.Close();
	dmEta = new TH1D("","",temp_eta->GetNbinsX(),temp_eta->GetXaxis()->GetXmin(),temp_eta->GetXaxis()->GetXmax());
	dmEta->Add(temp_eta);
	dmMolt = new TH1D("","",temp_mol->GetNbinsX(),temp_mol->GetXaxis()->GetXmin(),temp_mol->GetXaxis()->GetXmax());
	dmMolt->Add(temp_mol);
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
double MyRandom::RndTheta1() {
	double eta=dmEta->GetRandom();
	double temp=TMath::Exp(-eta);
	return 2.*TMath::ATan(temp);
}
