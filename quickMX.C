{
string pipORpim = "pim";

TFile *tf = new TFile("data.s1.n1500.BiSc9.MoCo11.__0000000000009000__.root");

int NxBins = 5;
int NQQBins = 2;
int NzBins = 4;
int NPT2Bins = 4;

TH1F *hMX[NxBins][NQQBins][NzBins][NPT2Bins];

TCanvas *can = new TCanvas("can", "can", 5, 5, 1000, 700);
can->Divide(NxBins, NQQBins, 0.00001, 0.00001);
//can->Divide(NxBins, NQQBins, 0.02, 0.02);

for(int x = 0; x < NxBins; x++) {
for(int QQ = 0; QQ < NQQBins; QQ++) {

//can->cd(NxBins*(NQQBins - QQ - 1) + x + 1)->Divide(NzBins, NPT2Bins, 0.00001, 0.00001);
can->cd(NxBins*(NQQBins - QQ - 1) + x + 1)->Divide(NzBins, NPT2Bins, 0, 0);

for(int z = 0; z < NzBins; z++) {
for(int PT2 = 0; PT2 < NPT2Bins; PT2++) {

hMX[x][QQ][z][PT2] = (TH1F*) tf->Get(Form("rec_%s_MX_x%i_QQ%i_z%i_PT2%i", pipORpim.c_str(), x, QQ, z, PT2));

can->cd(NxBins*(NQQBins - QQ - 1) + x + 1)->cd(NzBins*(NPT2Bins - PT2 - 1) + z + 1);

hMX[x][QQ][z][PT2]->Draw();

}}}}

}
