{
gStyle->SetOptStat(0);
TFile *sysSector = new TFile("/home/naharrison/mysidis-histos/Sector_systematics.root");

int nMBins = 0;
float MSysSum = 0.0;
int nAcBins = 0;
float AcSysSum = 0.0;
int nAccBins = 0;
float AccSysSum = 0.0;

string pipORpim = "pip";

for(int xBin = 0; xBin < 5; xBin++) {
for(int QQBin = 0; QQBin < 2; QQBin++) {
for(int zBin = 0; zBin < 18; zBin++) {
for(int PT2Bin = 0; PT2Bin < 20; PT2Bin++) {
if(!(xBin == 4 && QQBin == 1)) {

	TH1F *hSectorM = (TH1F*) sysSector->Get(Form("hM_%s_%i_%i_%i_%i", pipORpim.c_str(), xBin, QQBin, zBin, PT2Bin));
	TH1F *hSectorAc = (TH1F*) sysSector->Get(Form("hAc_%s_%i_%i_%i_%i", pipORpim.c_str(), xBin, QQBin, zBin, PT2Bin));
	TH1F *hSectorAcc = (TH1F*) sysSector->Get(Form("hAcc_%s_%i_%i_%i_%i", pipORpim.c_str(), xBin, QQBin, zBin, PT2Bin));
	
	// ### calculate the sector dependence systematic error contribution ###
	// ### from the 6 measurements (one for each sector) ###
	// ### Delta^2 = delta_stat^2 + delta_sys^2 ###
	float MAverage, AcAverage, AccAverage;
	MAverage = AcAverage = AccAverage = 0.0;
	for(int s = 0; s < 6; s++) {
		MAverage += hSectorM->GetBinContent(s+1);
		AcAverage += hSectorAc->GetBinContent(s+1);
		AccAverage += hSectorAcc->GetBinContent(s+1);
	}
	MAverage = MAverage/6.0;
	AcAverage = AcAverage/6.0;
	AccAverage = AccAverage/6.0;
	
	float MDelta, AcDelta, AccDelta;
	MDelta = AcDelta = AccDelta = 0.0;
	for(int s = 0; s < 6; s++) {
		MDelta += (hSectorM->GetBinContent(s+1) - MAverage)*(hSectorM->GetBinContent(s+1) - MAverage);
		AcDelta += (hSectorAc->GetBinContent(s+1) - AcAverage)*(hSectorAc->GetBinContent(s+1) - AcAverage);
		AccDelta += (hSectorAcc->GetBinContent(s+1) - AccAverage)*(hSectorAcc->GetBinContent(s+1) - AccAverage);
	}
	MDelta = sqrt(MDelta)/sqrt(6.0);
	AcDelta = sqrt(AcDelta)/sqrt(6.0);
	AccDelta = sqrt(AccDelta)/sqrt(6.0);
	
	float Mdelta_stat, Acdelta_stat, Accdelta_stat;
	Mdelta_stat = Acdelta_stat = Accdelta_stat = 0.0;
	for(int s = 0; s < 6; s++) {
		Mdelta_stat += hSectorM->GetBinError(s+1)*hSectorM->GetBinError(s+1);
		Acdelta_stat += hSectorAc->GetBinError(s+1)*hSectorAc->GetBinError(s+1);
		Accdelta_stat += hSectorAcc->GetBinError(s+1)*hSectorAcc->GetBinError(s+1);
	}
	Mdelta_stat = sqrt(Mdelta_stat)/sqrt(6.0);
	Acdelta_stat = sqrt(Acdelta_stat)/sqrt(6.0);
	Accdelta_stat = sqrt(Accdelta_stat)/sqrt(6.0);
	
	float Mdelta_sys, Acdelta_sys, Accdelta_sys;
	if(MDelta > Mdelta_stat) Mdelta_sys = sqrt(MDelta*MDelta - Mdelta_stat*Mdelta_stat);
	else Mdelta_sys = 0.0;
	if(AcDelta > Acdelta_stat) Acdelta_sys = sqrt(AcDelta*AcDelta - Acdelta_stat*Acdelta_stat);
	else Acdelta_sys = 0.0;
	if(AccDelta > Accdelta_stat) Accdelta_sys = sqrt(AccDelta*AccDelta - Accdelta_stat*Accdelta_stat);
	else Accdelta_sys = 0.0;
	
	// Get average sector systematic error
	// If stat error = 0 the bin is empty. Not counting these.
	if(Mdelta_stat > 0.0) {
		nMBins++;
		MSysSum = MSysSum + Mdelta_sys;
	}
	if(Acdelta_stat > 0.0) {
		nAcBins++;
		AcSysSum = AcSysSum + Acdelta_sys;
	}
	if(Accdelta_stat > 0.0) {
		nAccBins++;
		AccSysSum = AccSysSum + Accdelta_sys;
	}

	if(zBin == 0 && PT2Bin == 0) cout<<xBin<<" "<<QQBin<<endl;

}}}}}

cout<<endl;
cout<<MSysSum<<"/"<<nMBins<<endl;
cout<<AcSysSum<<"/"<<nAcBins<<endl;
cout<<AccSysSum<<"/"<<nAccBins<<endl;

}
