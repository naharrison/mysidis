{


int hapNumbers[4] = {0, 1, 2, 9};

TFile *rootFile = new TFile("MAcAcc.root", "recreate");

for(int ix = 0; ix < 5; ix++) {
for(int iQQ = 0; iQQ < 2; iQQ++) {
if(!(ix == 4 && iQQ == 1)) {
for(int iz = 0; iz < 18; iz++) {
for(int iPT2 = 0; iPT2 < 20; iPT2++) {
for(int iit = 0; iit < 3; iit++) {
for(int ihap : hapNumbers) {

	ifstream pipMfile(Form("/home/naharrison/MAcAccfiles/pip_M_it%i_hap%i_BiSc5_x%iQQ%iz%iPT2%i.txt", iit, ihap, ix, iQQ, iz, iPT2));
	if(pipMfile) {
		float pipM, pipMe;
		pipMfile>>pipM>>pipMe;
		TH1F *hpipM = new TH1F(Form("pip_M_it%i_hap%i_BiSc5_x%iQQ%iz%iPT2%i", iit, ihap, ix, iQQ, iz, iPT2), Form("pip_M_it%i_hap%i_BiSc5_x%iQQ%iz%iPT2%i", iit, ihap, ix, iQQ, iz, iPT2), 1, 0, 1);
		hpipM->SetBinContent(1, pipM);
		hpipM->SetBinError(1, pipMe);
	}
	pipMfile.close();

}}}}}}}

rootFile->Write();


}
