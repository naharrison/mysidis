{
TFile *MAcAccFile = new TFile("/home/naharrison/mysidis-histos/MAcAcc.root", "READ");
TFile *sysFile = new TFile("/home/naharrison/mysidis-histos/Total_systematics.root", "READ");

for(int xBin = 0; xBin < 5; xBin++) {
for(int QQBin = 0; QQBin < 2; QQBin++) {
for(int zBin = 0; zBin < 18; zBin++) {
for(int PT2Bin = 0; PT2Bin < 20; PT2Bin++) {
if(!(xBin == 4 && QQBin == 1)) {

	// pip:
	string pipSysHistName = Form("hMAcAccSys_pip_%i_%i_%i_%i", xBin, QQBin, zBin, PT2Bin);
	string pipMhistoName = Form("pip_M_it2_hap2_BiSc5_x%iQQ%iz%iPT2%i", xBin, QQBin, zBin, PT2Bin);
	string pipAcAcchistoName = Form("pip_AcAcc_it2_hap2_BiSc5_x%iQQ%iz%iPT2%i", xBin, QQBin, zBin, PT2Bin);
	
	if(sysFile->GetListOfKeys()->Contains(pipSysHistName.c_str())) {
	
		if(MAcAccFile->GetListOfKeys()->Contains(pipMhistoName.c_str())) {
			cout<<"pip_M_"<<xBin<<"_"<<QQBin<<"_"<<zBin<<"_"<<PT2Bin<<endl;
		}

		if(MAcAccFile->GetListOfKeys()->Contains(pipAcAcchistoName.c_str())) {
			cout<<"pip_AcAcc_"<<xBin<<"_"<<QQBin<<"_"<<zBin<<"_"<<PT2Bin<<endl;
		}
	
	}

	// pim:
	string pimSysHistName = Form("hMAcAccSys_pim_%i_%i_%i_%i", xBin, QQBin, zBin, PT2Bin);
	string pimMhistoName = Form("pim_M_it2_hap2_BiSc5_x%iQQ%iz%iPT2%i", xBin, QQBin, zBin, PT2Bin);
	string pimAcAcchistoName = Form("pim_AcAcc_it2_hap2_BiSc5_x%iQQ%iz%iPT2%i", xBin, QQBin, zBin, PT2Bin);
	
	if(sysFile->GetListOfKeys()->Contains(pimSysHistName.c_str())) {
	
		if(MAcAccFile->GetListOfKeys()->Contains(pimMhistoName.c_str())) {
			cout<<"pim_M_"<<xBin<<"_"<<QQBin<<"_"<<zBin<<"_"<<PT2Bin<<endl;
		}

		if(MAcAccFile->GetListOfKeys()->Contains(pimAcAcchistoName.c_str())) {
			cout<<"pim_AcAcc_"<<xBin<<"_"<<QQBin<<"_"<<zBin<<"_"<<PT2Bin<<endl;
		}
	
	}


}}}}}


}
