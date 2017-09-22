#include <fstream>
#include <iostream>
#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "programFiles/functions.C"
#include "programFiles/eID.C"
#include "programFiles/hadronID.C"
#include "programFiles/getGenIndices.C"
#include "MomCorr.C"

// for bin centering corrections... can also be use for other things
// inefficient code but good enough

void mysidis_BCC(int accIterationN = 2, int filestart = 1, int Nfiles = 10, int ExpOrSim = 0, bool do_momCorr_e = 1, bool do_momCorr_pions = 1, int e_zvertex_strict = 0, int e_ECsampling_strict = 0, int e_ECoutVin_strict = 0, int e_ECgeometric_strict = 0, int e_CCthetaMatching_strict = 0, int e_R1fid_strict = 0, int e_R3fid_strict = 0, int e_CCphiMatching_strict = 0, int e_CCfiducial_strict = 0, int yCut_strict = 0, int pip_vvp_strict = 0, int pip_R1fid_strict = 0, int pip_MXcut_strict = 0, int pim_vvp_strict = 0, int pim_R1fid_strict = 0, int pim_MXcut_strict = 0) // ExpOrSim = 0(MC) or 1(data) // note: momCorr is only done to data (even if it's turned on for MC)
{
TStopwatch *stopwat = new TStopwatch();

int MCversion = 12;

MomCorr_e1f *MomCorr = new MomCorr_e1f();

// %%%%% some cut values %%%%%
float WMin = 2.05;
float yMax = 0.85; // default value (for when yCut_strict = 0)
if(yCut_strict == 1) yMax = 0.8;
// %%%%% end some cut values %%%%%

// %%%%% x, QQ, z, PT2 bin definitions %%%%%
int binSchemeOpt = 55; // when you change this, there are two spots below that you have to update... 55 = bin scheme 5 shrunk into a small bin at the center of the original bin

int NxBins, NQQBins, NzBins_pip, NPT2Bins_pip, NzBins_pim, NPT2Bins_pim;
int NphihBins;
float phihMin, phihMax;

if(binSchemeOpt == 55)
{
NxBins = 5;
NQQBins = 2;
NzBins_pip = 18;
NPT2Bins_pip = 20;
NzBins_pim = 18;
NPT2Bins_pim = 20;
NphihBins = 36*5; // easier to just multiply this one by 5
phihMin = -180;
phihMax = 180;
}

//float xLimits[NxBins + 1], QQLimits[NxBins][NQQBins + 1], zLimits_pip[NzBins_pip + 1], PT2Limits_pip[NPT2Bins_pip + 1], zLimits_pim[NzBins_pim + 1], PT2Limits_pim[NPT2Bins_pim + 1];
float xLow[NxBins], xHigh[NxBins], QQLow[NxBins][NQQBins], QQHigh[NxBins][NQQBins], zLow_pip[NzBins_pip], zHigh_pip[NzBins_pip], PT2Low_pip[NzBins_pip], PT2High_pip[NzBins_pip], zLow_pim[NzBins_pim], zHigh_pim[NzBins_pim], PT2Low_pim[NzBins_pim], PT2High_pim[NzBins_pim];

if(binSchemeOpt == 55)
{
xLow[0] = 0.17 - 0.01;
xLow[1] = 0.25 - 0.01;
xLow[2] = 0.35 - 0.01;
xLow[3] = 0.42 - 0.01;
xLow[4] = 0.53 - 0.01;

xHigh[0] = 0.17 + 0.01;
xHigh[1] = 0.25 + 0.01;
xHigh[2] = 0.35 + 0.01;
xHigh[3] = 0.42 + 0.01;
xHigh[4] = 0.53 + 0.01;

QQLow[0][0] = 1.15 - 0.05;
QQLow[0][1] = 1.35 - 0.05;
QQLow[1][0] = 1.45 - 0.05;
QQLow[1][1] = 1.95 - 0.05;
QQLow[2][0] = 2.05 - 0.05;
QQLow[2][1] = 2.55 - 0.05;
QQLow[3][0] = 2.65 - 0.05;
QQLow[3][1] = 3.15 - 0.05;
QQLow[4][0] = 4.05 - 0.05;
QQLow[4][1] = 5.55 - 0.05;

QQHigh[0][0] = 1.15 + 0.05;
QQHigh[0][1] = 1.35 + 0.05;
QQHigh[1][0] = 1.45 + 0.05;
QQHigh[1][1] = 1.95 + 0.05;
QQHigh[2][0] = 2.05 + 0.05;
QQHigh[2][1] = 2.55 + 0.05;
QQHigh[3][0] = 2.65 + 0.05;
QQHigh[3][1] = 3.15 + 0.05;
QQHigh[4][0] = 4.05 + 0.05;
QQHigh[4][1] = 5.55 + 0.05;

zLow_pip[0] = 0.025 - 0.005;
zLow_pip[1] = 0.075 - 0.005;
zLow_pip[2] = 0.125 - 0.005;
zLow_pip[3] = 0.175 - 0.005;
zLow_pip[4] = 0.225 - 0.005;
zLow_pip[5] = 0.275 - 0.005;
zLow_pip[6] = 0.325 - 0.005;
zLow_pip[7] = 0.375 - 0.005;
zLow_pip[8] = 0.425 - 0.005;
zLow_pip[9] = 0.475 - 0.005;
zLow_pip[10] = 0.525 - 0.005;
zLow_pip[11] = 0.575 - 0.005;
zLow_pip[12] = 0.625 - 0.005;
zLow_pip[13] = 0.675 - 0.005;
zLow_pip[14] = 0.725 - 0.005;
zLow_pip[15] = 0.775 - 0.005;
zLow_pip[16] = 0.825 - 0.005;
zLow_pip[17] = 0.875 - 0.005;

zHigh_pip[0] = 0.025 + 0.005;
zHigh_pip[1] = 0.075 + 0.005;
zHigh_pip[2] = 0.125 + 0.005;
zHigh_pip[3] = 0.175 + 0.005;
zHigh_pip[4] = 0.225 + 0.005;
zHigh_pip[5] = 0.275 + 0.005;
zHigh_pip[6] = 0.325 + 0.005;
zHigh_pip[7] = 0.375 + 0.005;
zHigh_pip[8] = 0.425 + 0.005;
zHigh_pip[9] = 0.475 + 0.005;
zHigh_pip[10] = 0.525 + 0.005;
zHigh_pip[11] = 0.575 + 0.005;
zHigh_pip[12] = 0.625 + 0.005;
zHigh_pip[13] = 0.675 + 0.005;
zHigh_pip[14] = 0.725 + 0.005;
zHigh_pip[15] = 0.775 + 0.005;
zHigh_pip[16] = 0.825 + 0.005;
zHigh_pip[17] = 0.875 + 0.005;

PT2Low_pip[0] = 0.025 - 0.005;
PT2Low_pip[1] = 0.075 - 0.005;
PT2Low_pip[2] = 0.125 - 0.005;
PT2Low_pip[3] = 0.175 - 0.005;
PT2Low_pip[4] = 0.225 - 0.005;
PT2Low_pip[5] = 0.275 - 0.005;
PT2Low_pip[6] = 0.325 - 0.005;
PT2Low_pip[7] = 0.375 - 0.005;
PT2Low_pip[8] = 0.425 - 0.005;
PT2Low_pip[9] = 0.475 - 0.005;
PT2Low_pip[10] = 0.525 - 0.005;
PT2Low_pip[11] = 0.575 - 0.005;
PT2Low_pip[12] = 0.625 - 0.005;
PT2Low_pip[13] = 0.675 - 0.005;
PT2Low_pip[14] = 0.725 - 0.005;
PT2Low_pip[15] = 0.775 - 0.005;
PT2Low_pip[16] = 0.825 - 0.005;
PT2Low_pip[17] = 0.875 - 0.005;
PT2Low_pip[18] = 0.925 - 0.005;
PT2Low_pip[19] = 0.975 - 0.005;

PT2High_pip[0] = 0.025 + 0.005;
PT2High_pip[1] = 0.075 + 0.005;
PT2High_pip[2] = 0.125 + 0.005;
PT2High_pip[3] = 0.175 + 0.005;
PT2High_pip[4] = 0.225 + 0.005;
PT2High_pip[5] = 0.275 + 0.005;
PT2High_pip[6] = 0.325 + 0.005;
PT2High_pip[7] = 0.375 + 0.005;
PT2High_pip[8] = 0.425 + 0.005;
PT2High_pip[9] = 0.475 + 0.005;
PT2High_pip[10] = 0.525 + 0.005;
PT2High_pip[11] = 0.575 + 0.005;
PT2High_pip[12] = 0.625 + 0.005;
PT2High_pip[13] = 0.675 + 0.005;
PT2High_pip[14] = 0.725 + 0.005;
PT2High_pip[15] = 0.775 + 0.005;
PT2High_pip[16] = 0.825 + 0.005;
PT2High_pip[17] = 0.875 + 0.005;
PT2High_pip[18] = 0.925 + 0.005;
PT2High_pip[19] = 0.975 + 0.005;

zLow_pim[0] = 0.025 - 0.005;
zLow_pim[1] = 0.075 - 0.005;
zLow_pim[2] = 0.125 - 0.005;
zLow_pim[3] = 0.175 - 0.005;
zLow_pim[4] = 0.225 - 0.005;
zLow_pim[5] = 0.275 - 0.005;
zLow_pim[6] = 0.325 - 0.005;
zLow_pim[7] = 0.375 - 0.005;
zLow_pim[8] = 0.425 - 0.005;
zLow_pim[9] = 0.475 - 0.005;
zLow_pim[10] = 0.525 - 0.005;
zLow_pim[11] = 0.575 - 0.005;
zLow_pim[12] = 0.625 - 0.005;
zLow_pim[13] = 0.675 - 0.005;
zLow_pim[14] = 0.725 - 0.005;
zLow_pim[15] = 0.775 - 0.005;
zLow_pim[16] = 0.825 - 0.005;
zLow_pim[17] = 0.875 - 0.005;

zHigh_pim[0] = 0.025 + 0.005;
zHigh_pim[1] = 0.075 + 0.005;
zHigh_pim[2] = 0.125 + 0.005;
zHigh_pim[3] = 0.175 + 0.005;
zHigh_pim[4] = 0.225 + 0.005;
zHigh_pim[5] = 0.275 + 0.005;
zHigh_pim[6] = 0.325 + 0.005;
zHigh_pim[7] = 0.375 + 0.005;
zHigh_pim[8] = 0.425 + 0.005;
zHigh_pim[9] = 0.475 + 0.005;
zHigh_pim[10] = 0.525 + 0.005;
zHigh_pim[11] = 0.575 + 0.005;
zHigh_pim[12] = 0.625 + 0.005;
zHigh_pim[13] = 0.675 + 0.005;
zHigh_pim[14] = 0.725 + 0.005;
zHigh_pim[15] = 0.775 + 0.005;
zHigh_pim[16] = 0.825 + 0.005;
zHigh_pim[17] = 0.875 + 0.005;

PT2Low_pim[0] = 0.025 - 0.005;
PT2Low_pim[1] = 0.075 - 0.005;
PT2Low_pim[2] = 0.125 - 0.005;
PT2Low_pim[3] = 0.175 - 0.005;
PT2Low_pim[4] = 0.225 - 0.005;
PT2Low_pim[5] = 0.275 - 0.005;
PT2Low_pim[6] = 0.325 - 0.005;
PT2Low_pim[7] = 0.375 - 0.005;
PT2Low_pim[8] = 0.425 - 0.005;
PT2Low_pim[9] = 0.475 - 0.005;
PT2Low_pim[10] = 0.525 - 0.005;
PT2Low_pim[11] = 0.575 - 0.005;
PT2Low_pim[12] = 0.625 - 0.005;
PT2Low_pim[13] = 0.675 - 0.005;
PT2Low_pim[14] = 0.725 - 0.005;
PT2Low_pim[15] = 0.775 - 0.005;
PT2Low_pim[16] = 0.825 - 0.005;
PT2Low_pim[17] = 0.875 - 0.005;
PT2Low_pim[18] = 0.925 - 0.005;
PT2Low_pim[19] = 0.975 - 0.005;

PT2High_pim[0] = 0.025 + 0.005;
PT2High_pim[1] = 0.075 + 0.005;
PT2High_pim[2] = 0.125 + 0.005;
PT2High_pim[3] = 0.175 + 0.005;
PT2High_pim[4] = 0.225 + 0.005;
PT2High_pim[5] = 0.275 + 0.005;
PT2High_pim[6] = 0.325 + 0.005;
PT2High_pim[7] = 0.375 + 0.005;
PT2High_pim[8] = 0.425 + 0.005;
PT2High_pim[9] = 0.475 + 0.005;
PT2High_pim[10] = 0.525 + 0.005;
PT2High_pim[11] = 0.575 + 0.005;
PT2High_pim[12] = 0.625 + 0.005;
PT2High_pim[13] = 0.675 + 0.005;
PT2High_pim[14] = 0.725 + 0.005;
PT2High_pim[15] = 0.775 + 0.005;
PT2High_pim[16] = 0.825 + 0.005;
PT2High_pim[17] = 0.875 + 0.005;
PT2High_pim[18] = 0.925 + 0.005;
PT2High_pim[19] = 0.975 + 0.005;
}
// %%%%% end x, QQ, z, PT2 bin definitions %%%%%

// %%%%% read the input files into a TChain %%%%%
int NtotalFiles = 11625;
if(ExpOrSim == 0 && MCversion == 3) NtotalFiles = 33471;
if(ExpOrSim == 0 && MCversion == 8) NtotalFiles = 32171;
if(ExpOrSim == 0 && MCversion == 9) NtotalFiles = 1995;
if(ExpOrSim == 0 && MCversion == 10) NtotalFiles = 3950;
if(ExpOrSim == 0 && MCversion == 11) NtotalFiles = 3931;
if(ExpOrSim == 0 && MCversion == 12) NtotalFiles = 32255;

string firstfilename = "";
string lastfilename = "";

ifstream filelist;
if(ExpOrSim == 0) filelist.open(Form("programFiles/v%i_MCFiles.txt", MCversion));
else filelist.open("programFiles/dataFiles.txt");

TChain *h22chain = new TChain("h22");

int kStop = Nfiles + filestart;
if(kStop > NtotalFiles+1) kStop = NtotalFiles+1;
for(int k = 1; k < kStop; k++)
{
string filename;
filelist>>filename;
if(k >= filestart) h22chain->Add(filename.c_str());
if(k == filestart) firstfilename = filename;
if(k == kStop - 1) lastfilename = filename;
}

if(ExpOrSim == 1)
{
firstfilename = firstfilename.substr(36,9); // trim the string down so it's just the run#.subrun# (e.g. 38458.a00)
lastfilename = lastfilename.substr(36,9);
}
cout<<"first file: "<<firstfilename<<" ... last file: "<<lastfilename<<endl;
// %%%%% end read the input files into a TChain %%%%%

// %%%%% define variables %%%%%
Float_t e_mass = 0.000511; // GeV
Float_t prot_mass = 0.938272; // GeV
Float_t pip_mass = 0.13957; // GeV
Float_t pim_mass = 0.13957; // GeV
Float_t speed_of_light = 29.9792458; // cm/ns
Float_t Beam_Energy = 5.498; // GeV
Float_t pi = 3.14159265359;
Float_t pi180 = pi/180.0;
Float_t pi180_inv = 180.0/pi;

Int_t mcnentr;
Int_t mcid[35];
Float_t mctheta[35];
Float_t mcphi[35];
Float_t mcp[35];

Int_t gpart;

Float_t p[35];
Int_t q[35];
Float_t cx[35];
Float_t cy[35];
Float_t cz[35];
Float_t vz[35];
Float_t vy[35];
Float_t vx[35];
Float_t b[35];

UChar_t dc_sect[35];
Float_t tl1_x[35];
Float_t tl1_y[35];
Float_t tl3_x[35]; // used to be dc_xsc in h10
Float_t tl3_y[35];
Float_t tl3_z[35];
Float_t tl3_cx[35]; // used to be dc_cxsc in h10
Float_t tl3_cy[35];
Float_t tl3_cz[35];

UChar_t ec_sect[35];
Float_t etot[35];
Float_t ec_ei[35];
Float_t ec_eo[35];
Float_t ec_t[35];
Float_t ech_x[35];
Float_t ech_y[35];
Float_t ech_z[35];

UChar_t sc_sect[35];
Float_t sc_t[35];
Float_t sc_r[35];
UChar_t sc_pd[35];

UChar_t cc_sect[35];
UShort_t cc_segm[35];
UShort_t nphe[35];
// %%%%% end define variables %%%%%

// %%%%% set branch addresses %%%%%
if(ExpOrSim == 0)
{
h22chain->SetBranchAddress("mcnentr", &mcnentr);
h22chain->SetBranchAddress("mcid", mcid);
h22chain->SetBranchAddress("mctheta", mctheta);
h22chain->SetBranchAddress("mcphi", mcphi);
h22chain->SetBranchAddress("mcp", mcp);
}

h22chain->SetBranchAddress("gpart", &gpart);
h22chain->SetBranchAddress("p", p);
h22chain->SetBranchAddress("q", q);
h22chain->SetBranchAddress("cx", cx);
h22chain->SetBranchAddress("cy", cy);
h22chain->SetBranchAddress("cz", cz);
h22chain->SetBranchAddress("vz", vz);
h22chain->SetBranchAddress("vy", vy);
h22chain->SetBranchAddress("vx", vx);
h22chain->SetBranchAddress("b", b);

h22chain->SetBranchAddress("dc_sect", dc_sect);
h22chain->SetBranchAddress("tl1_x", tl1_x);
h22chain->SetBranchAddress("tl1_y", tl1_y);
h22chain->SetBranchAddress("tl3_x", tl3_x);
h22chain->SetBranchAddress("tl3_y", tl3_y);
h22chain->SetBranchAddress("tl3_z", tl3_z);
h22chain->SetBranchAddress("tl3_cx", tl3_cx);
h22chain->SetBranchAddress("tl3_cy", tl3_cy);
h22chain->SetBranchAddress("tl3_cz", tl3_cz);

h22chain->SetBranchAddress("ec_sect", ec_sect);
h22chain->SetBranchAddress("etot", etot);
h22chain->SetBranchAddress("ec_ei", ec_ei);
h22chain->SetBranchAddress("ec_eo", ec_eo);
h22chain->SetBranchAddress("ec_t", ec_t);
h22chain->SetBranchAddress("ech_x", ech_x);
h22chain->SetBranchAddress("ech_y", ech_y);
h22chain->SetBranchAddress("ech_z", ech_z);

h22chain->SetBranchAddress("sc_sect", sc_sect);
h22chain->SetBranchAddress("sc_t", sc_t);
h22chain->SetBranchAddress("sc_r", sc_r);
h22chain->SetBranchAddress("sc_pd", sc_pd);

h22chain->SetBranchAddress("cc_sect", cc_sect);
h22chain->SetBranchAddress("cc_segm", cc_segm);
h22chain->SetBranchAddress("nphe", nphe);
// %%%%% end set branch addresses %%%%%

// %%%%% stuff for output root files %%%%%
string outfilename;
TFile *outputfile;
if(ExpOrSim == 1) outfilename = Form("data_BCC.s%i.n%i.BiSc%i.MoCo%i%i.__%i%i%i%i%i%i%i%i%i%i%i%i%i%i%i%i__.root", filestart, Nfiles, binSchemeOpt, do_momCorr_e, do_momCorr_pions, e_zvertex_strict, e_ECsampling_strict, e_ECoutVin_strict, e_ECgeometric_strict, e_CCthetaMatching_strict, e_R1fid_strict, e_R3fid_strict, e_CCphiMatching_strict, e_CCfiducial_strict, yCut_strict, pip_vvp_strict, pip_R1fid_strict, pip_MXcut_strict, pim_vvp_strict, pim_R1fid_strict, pim_MXcut_strict);
if(ExpOrSim == 0) outfilename = Form("MonteCarlo_v%i_BCC.it%i.s%i.n%i.BiSc%i.__%i%i%i%i%i%i%i%i%i%i%i%i%i%i%i%i__.root", MCversion, accIterationN, filestart, Nfiles, binSchemeOpt, e_zvertex_strict, e_ECsampling_strict, e_ECoutVin_strict, e_ECgeometric_strict, e_CCthetaMatching_strict, e_R1fid_strict, e_R3fid_strict, e_CCphiMatching_strict, e_CCfiducial_strict, yCut_strict, pip_vvp_strict, pip_R1fid_strict, pip_MXcut_strict, pim_vvp_strict, pim_R1fid_strict, pim_MXcut_strict);
cout<<"output filename is "<<outfilename<<endl;
outputfile = new TFile(outfilename.c_str(),"recreate");
// %%%%% end stuff for output root files %%%%%

// %%%%% define histograms %%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TH1F *h_pip_phih[NxBins][NQQBins][NzBins_pip][NPT2Bins_pip][2]; // 2 for gen(0) or rec(1)

for(int x = 0; x < NxBins; x++){
for(int QQ = 0; QQ < NQQBins; QQ++){
for(int z = 0; z < NzBins_pip; z++){
for(int PT2 = 0; PT2 < NPT2Bins_pip; PT2++){
if(ExpOrSim == 0)
{
h_pip_phih[x][QQ][z][PT2][0] = new TH1F(Form("gen_pip_phih_x%i_QQ%i_z%i_PT2%i",x,QQ,z,PT2), Form("gen_pip_phih_x%i_QQ%i_z%i_PT2%i",x,QQ,z,PT2), NphihBins, phihMin, phihMax);
}
h_pip_phih[x][QQ][z][PT2][1] = new TH1F(Form("rec_pip_phih_x%i_QQ%i_z%i_PT2%i",x,QQ,z,PT2), Form("rec_pip_phih_x%i_QQ%i_z%i_PT2%i",x,QQ,z,PT2), NphihBins, phihMin, phihMax);
}}}}

TH1F *h_pim_phih[NxBins][NQQBins][NzBins_pim][NPT2Bins_pim][2]; // 2 for gen(0) or rec(1)

for(int x = 0; x < NxBins; x++){
for(int QQ = 0; QQ < NQQBins; QQ++){
for(int z = 0; z < NzBins_pim; z++){
for(int PT2 = 0; PT2 < NPT2Bins_pim; PT2++){
if(ExpOrSim == 0)
{
h_pim_phih[x][QQ][z][PT2][0] = new TH1F(Form("gen_pim_phih_x%i_QQ%i_z%i_PT2%i",x,QQ,z,PT2), Form("gen_pim_phih_x%i_QQ%i_z%i_PT2%i",x,QQ,z,PT2), NphihBins, phihMin, phihMax);
}
h_pim_phih[x][QQ][z][PT2][1] = new TH1F(Form("rec_pim_phih_x%i_QQ%i_z%i_PT2%i",x,QQ,z,PT2), Form("rec_pim_phih_x%i_QQ%i_z%i_PT2%i",x,QQ,z,PT2), NphihBins, phihMin, phihMax);
}}}}
// %%%%% end define histograms %%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// %%%%% POTATOES %%%%%
// %%%%%%%%%%%%%%%%%%%%
cout<<"entries: "<<h22chain->GetEntries()/1000<<" thousand"<<endl;

for(int i = 0; i < h22chain->GetEntries(); i++) // loop over entries
{
h22chain->GetEntry(i);

bool genExclusiveEvent = 0; // if the generated event is an exclusive event, i.e. ep --> e pi+ n, we don't want to analyze it at all (the generator give unrealistic kinematics)
if(ExpOrSim == 0 && mcnentr == 3)
{
if((mcid[1] == 211 && mcid[2] == 2112) || (mcid[1] == 2112 && mcid[2] == 211)) genExclusiveEvent = 1; // the first (zeroth) generated particle is always the scattered electron
}

if(genExclusiveEvent == 0)
{

string currentrunno_string = "-123";
if(ExpOrSim == 1)
{
currentrunno_string = h22chain->GetCurrentFile()->GetName();
currentrunno_string = currentrunno_string.substr(36,5);
}
int currentrunno = atoi(currentrunno_string.c_str()); // converts the string to an int

TLorentzVector V4k(0.0, 0.0, Beam_Energy, Beam_Energy); // x, y, z, t
TLorentzVector V4ISproton(0.0, 0.0, 0.0, prot_mass); // IS = Initial State
int e_index[2] = {-123, -123};
int pip_index[2] = {-123, -123};
int pim_index[2] = {-123, -123};
int prot_index[2] = {-123, -123};

// %%%%% get gen. particle indices  %%%%%
if(ExpOrSim == 0)
{
e_index[0] = 0; // the first (zeroth) generated particle is always the scattered electron
vector<int> genIndices = getGenIndices(mcnentr, mcid, mcp);
pip_index[0] = genIndices[0];
pim_index[0] = genIndices[1];
prot_index[0] = genIndices[2];
}
// %%%%% end get gen. particle indices  %%%%%

// %%%%% electron ID %%%%%
e_index[1] = eID(gpart, q, p, cc_sect, sc_sect, ec_sect, dc_sect, cx, cy, cz, tl1_x, tl1_y, tl3_x, tl3_y, tl3_z, tl3_cx, tl3_cy, tl3_cz, e_zvertex_strict, vz, vy, vx, e_ECsampling_strict, ExpOrSim, etot, e_ECoutVin_strict, ec_ei, ech_x, ech_y, ech_z, e_CCthetaMatching_strict, cc_segm, e_ECgeometric_strict, e_R1fid_strict, e_R3fid_strict, e_CCphiMatching_strict, sc_pd, e_CCfiducial_strict);
// %%%%% end electron ID %%%%%

// %%%%% set e- 3 and 4 vectors %%%%%
// can't set the vectors in a loop because gen uses p, theta, phi but rec uses p, cx, cy, cz
TVector3 V3_e[2]; // 2 for gen(0) and rec(1)
TLorentzVector V4_e[2];

if(e_index[0] > -122)
{
V3_e[0].SetXYZ(mcp[e_index[0]]*cos(pi180*mcphi[e_index[0]])*sin(pi180*mctheta[e_index[0]]), mcp[e_index[0]]*sin(pi180*mcphi[e_index[0]])*sin(pi180*mctheta[e_index[0]]), mcp[e_index[0]]*cos(pi180*mctheta[e_index[0]]));
V4_e[0].SetXYZT(V3_e[0].X(), V3_e[0].Y(), V3_e[0].Z(), sqrt(V3_e[0].Mag2() + pow(e_mass,2)));
}

if(e_index[1] > -122)
{
V3_e[1].SetXYZ(p[e_index[1]]*cx[e_index[1]], p[e_index[1]]*cy[e_index[1]], p[e_index[1]]*cz[e_index[1]]);
V4_e[1].SetXYZT(V3_e[1].X(), V3_e[1].Y(), V3_e[1].Z(), sqrt(V3_e[1].Mag2() + pow(e_mass,2)));
if(do_momCorr_e && ExpOrSim == 1) V4_e[1] = MomCorr->PcorN(V4_e[1], -1, 11);
}

TVector3 V3_H[2];
TLorentzVector V4_q[2], V4_H[2];

for(int j = 1; j >= ExpOrSim; j--)
{
V4_q[j] = V4k - V4_e[j];
V4_H[j] = V4_q[j] + V4ISproton;
V3_H[j].SetXYZ(V4_H[j].X(), V4_H[j].Y(), V4_H[j].Z());
}
// %%%%% end set e- 3 and 4 vectors %%%%%

// %%%%% hadron ID %%%%%
if(e_index[1] >= 0)
{
vector<int> hadronIndices = hadronID(gpart, e_index, q, p, sc_sect, dc_sect, sc_t, sc_r, sc_pd, pip_vvp_strict, pip_R1fid_strict, pip_MXcut_strict, ExpOrSim, ec_ei, ec_sect, cc_sect, nphe, ec_eo, cx, cy, cz, b, tl1_x, tl1_y, mcp, mcphi, mctheta, V4_H, currentrunno, pim_vvp_strict, pim_R1fid_strict, pim_MXcut_strict);
pip_index[1] = hadronIndices[0];
pim_index[1] = hadronIndices[1];
prot_index[1] = hadronIndices[2];
}
// %%%%% end hadron ID %%%%%

// %%%%%%%%% set remaining 4-vectors (and 3-vectors)  %%%%%%%%%%%
// can't set the vectors in a loop because gen uses p, theta, phi but rec uses p, cx, cy, cz
TVector3 V3_pip[2], V3_pim[2], V3_prot[2]; // 2 for gen(0) and rec(1)
TLorentzVector V4_pip[2], V4_pim[2], V4_prot[2];

if(pip_index[0] > -122)
{
V3_pip[0].SetXYZ(mcp[pip_index[0]]*cos(pi180*mcphi[pip_index[0]])*sin(pi180*mctheta[pip_index[0]]), mcp[pip_index[0]]*sin(pi180*mcphi[pip_index[0]])*sin(pi180*mctheta[pip_index[0]]), mcp[pip_index[0]]*cos(pi180*mctheta[pip_index[0]]));
V4_pip[0].SetXYZT(V3_pip[0].X(), V3_pip[0].Y(), V3_pip[0].Z(), sqrt(V3_pip[0].Mag2() + pow(pip_mass,2)));
}
if(pim_index[0] > -122)
{
V3_pim[0].SetXYZ(mcp[pim_index[0]]*cos(pi180*mcphi[pim_index[0]])*sin(pi180*mctheta[pim_index[0]]), mcp[pim_index[0]]*sin(pi180*mcphi[pim_index[0]])*sin(pi180*mctheta[pim_index[0]]), mcp[pim_index[0]]*cos(pi180*mctheta[pim_index[0]]));
V4_pim[0].SetXYZT(V3_pim[0].X(), V3_pim[0].Y(), V3_pim[0].Z(), sqrt(V3_pim[0].Mag2() + pow(pim_mass,2)));
}
if(prot_index[0] > -122)
{
V3_prot[0].SetXYZ(mcp[prot_index[0]]*cos(pi180*mcphi[prot_index[0]])*sin(pi180*mctheta[prot_index[0]]), mcp[prot_index[0]]*sin(pi180*mcphi[prot_index[0]])*sin(pi180*mctheta[prot_index[0]]), mcp[prot_index[0]]*cos(pi180*mctheta[prot_index[0]]));
V4_prot[0].SetXYZT(V3_prot[0].X(), V3_prot[0].Y(), V3_prot[0].Z(), sqrt(V3_prot[0].Mag2() + pow(prot_mass,2)));
}

if(pip_index[1] > -122)
{
V3_pip[1].SetXYZ(p[pip_index[1]]*cx[pip_index[1]], p[pip_index[1]]*cy[pip_index[1]], p[pip_index[1]]*cz[pip_index[1]]);
V4_pip[1].SetXYZT(V3_pip[1].X(), V3_pip[1].Y(), V3_pip[1].Z(), sqrt(V3_pip[1].Mag2() + pow(pip_mass,2)));
if(do_momCorr_pions && ExpOrSim == 1) V4_pip[1] = MomCorr->PcorN(V4_pip[1], 1, 211);
}
if(pim_index[1] > -122)
{
V3_pim[1].SetXYZ(p[pim_index[1]]*cx[pim_index[1]], p[pim_index[1]]*cy[pim_index[1]], p[pim_index[1]]*cz[pim_index[1]]);
V4_pim[1].SetXYZT(V3_pim[1].X(), V3_pim[1].Y(), V3_pim[1].Z(), sqrt(V3_pim[1].Mag2() + pow(pim_mass,2)));
if(do_momCorr_pions && ExpOrSim == 1) V4_pim[1] = MomCorr->PcorN(V4_pim[1], -1, -211);
}
if(prot_index[1] > -122)
{
V3_prot[1].SetXYZ(p[prot_index[1]]*cx[prot_index[1]], p[prot_index[1]]*cy[prot_index[1]], p[prot_index[1]]*cz[prot_index[1]]);
V4_prot[1].SetXYZT(V3_prot[1].X(), V3_prot[1].Y(), V3_prot[1].Z(), sqrt(V3_prot[1].Mag2() + pow(prot_mass,2)));
}

TLorentzVector V4_pip_prime[2], V4_pim_prime[2], V4_prot_prime[2];
TLorentzVector V4_X_epipX[2];
TLorentzVector V4_X_epimX[2];

for(int j = 1; j >= ExpOrSim; j--)
{
V4_q[j] = V4k - V4_e[j];
V4_H[j] = V4_q[j] + V4ISproton;
V3_H[j].SetXYZ(V4_H[j].X(), V4_H[j].Y(), V4_H[j].Z());

V4_pip_prime[j] = V4_pip[j];
V4_pip_prime[j].RotateZ(-1.0*V4_q[j].Phi() - pi);
V4_pip_prime[j].RotateY(V4_q[j].Theta());
V4_pip_prime[j].Boost(0.0, 0.0, V3_H[j].Mag()/V4_H[j].T());

V4_pim_prime[j] = V4_pim[j];
V4_pim_prime[j].RotateZ(-1.0*V4_q[j].Phi() - pi);
V4_pim_prime[j].RotateY(V4_q[j].Theta());
V4_pim_prime[j].Boost(0.0, 0.0, V3_H[j].Mag()/V4_H[j].T());

V4_prot_prime[j] = V4_prot[j];
V4_prot_prime[j].RotateZ(-1.0*V4_q[j].Phi() - pi);
V4_prot_prime[j].RotateY(V4_q[j].Theta());
V4_prot_prime[j].Boost(0.0, 0.0, V3_H[j].Mag()/V4_H[j].T());

V4_X_epipX[j] = V4_H[j] - V4_pip[j];
V4_X_epimX[j] = V4_H[j] - V4_pim[j];
}
// %%%%%%%%% end set remaining 4-vectors (and 3-vectors)  %%%%%%%%%%%

// %%%%%%%%% analysis of kinematics %%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
float W[2], y[2], QQ[2], x[2], pip_z[2], pip_PT2[2], pip_phih[2], MX_epipX[2];
int xBin[2] = {-1234, -1234};
int QQBin[2] = {-1234, -1234};
int pip_zBin[2] = {-1234, -1234};
int pip_PT2Bin[2] = {-1234, -1234};

float pim_z[2], pim_PT2[2], pim_phih[2], MX_epimX[2];
int pim_zBin[2] = {-1234, -1234};
int pim_PT2Bin[2] = {-1234, -1234};

// set some values:
for(int j = 1; j >= ExpOrSim; j--)
{
W[j] = V4_H[j].Mag();
y[j] = V4_q[j].T()/Beam_Energy;
QQ[j] = -V4_q[j]*V4_q[j];
x[j] = QQ[j]/(2.0*V4_q[j].T()*prot_mass);

pip_z[j] = V4_pip[j].T()/V4_q[j].T();
pip_PT2[j] = V4_pip_prime[j].Perp()*V4_pip_prime[j].Perp();
pip_phih[j] = pi180_inv*V4_pip_prime[j].Phi();
MX_epipX[j] = V4_X_epipX[j].Mag();

pim_z[j] = V4_pim[j].T()/V4_q[j].T();
pim_PT2[j] = V4_pim_prime[j].Perp()*V4_pim_prime[j].Perp();
pim_phih[j] = pi180_inv*V4_pim_prime[j].Phi();
MX_epimX[j] = V4_X_epimX[j].Mag();

if(pip_index[j] >= 0 || pim_index[j] >= 0)
{
xBin[j] = getBinN3(x[j], NxBins, xLow, xHigh);
QQBin[j] = getBinN3(QQ[j], NQQBins, QQLow[xBin[j]], QQHigh[xBin[j]]);

if(xBin[j] > -1 && QQBin[j] > -1)
{
if(pip_index[j] >= 0)
{
pip_zBin[j] = getBinN3(pip_z[j], NzBins_pip, zLow_pip, zHigh_pip);
pip_PT2Bin[j] = getBinN3(pip_PT2[j], NPT2Bins_pip, PT2Low_pip, PT2High_pip);
}

if(pim_index[j] >= 0)
{
pim_zBin[j] = getBinN3(pim_z[j], NzBins_pim, zLow_pim, zHigh_pim);
pim_PT2Bin[j] = getBinN3(pim_PT2[j], NPT2Bins_pim, PT2Low_pim, PT2High_pim);
}
}
}
}

/////// get weight factors: ///////
///////////////////////////////////
float pip_weight = 1.0;
float pim_weight = 1.0;

// iteration 1 = default model in haprad:
// pip:
if(ExpOrSim == 0 && accIterationN == 1 && xBin[0] > -1 && QQBin[0] > -1 && W[0] > WMin && y[0] < yMax && pip_index[0] >= 0)
{
ifstream term0file(Form("/scratch/dflt_term0_fromHAPRAD/dflt_term0_fromHAPRAD_BiSc%i_x%iQQ%iz%iPT2%i.txt", binSchemeOpt, xBin[0], QQBin[0], pip_zBin[0], pip_PT2Bin[0]));
if(term0file)
{
float term0;
term0file>>term0;

// from semi_inclusive_model.f: "Parton transverse momentum distribution (Gaussian fit)"
float ac = 1.2025E-10;
float bc = -5.2703E-02;
float cc = 3.7467E-01;
float dc = 6.5397E-02;
float ec = -2.2136E-01;
float fc = -1.0621E-01;

float sgmpt = ac + bc*x[0] + cc*pip_z[0] + dc*x[0]*x[0] + ec*pip_z[0]*pip_z[0] + fc*x[0]*pip_z[0];
if(sgmpt < 0.02) sgmpt = 0.02;
if(sgmpt > 0.15) sgmpt = 0.15;
float GTMD = exp((-(pip_PT2[0] + 0.0))/(2.0*sgmpt))*(1.0/(2.0*pi*sgmpt)); // the 0.0 part comes from "2*pl^2" in haprad, but "pl" seems to always be zero
float pi_thresh = sqrt(1.0 - ((prot_mass + pip_mass)*(prot_mass + pip_mass))/pow(MX_epipX[0], 2.0));
float rt = 1.0 - y[0] - (prot_mass*x[0]*y[0])/(2.0*Beam_Energy);
float rtz = sqrt(rt/(1.0 + ((2.0*prot_mass*x[0])/(y[0]*Beam_Energy))));
float h3 = -0.36544E-03*pow(x[0], -2.1855)*pow(1.0 - x[0], 3.4176)*pow(pip_z[0], -1.7567)*pow(1.0 - pip_z[0], 1.1272)*pow(log(QQ[0]/(0.25*0.25))/log(1.0/(0.25*0.25)), 8.9985);
float h4 = 0.10908E-02*pow(x[0], -0.35265E-06)*pow(1.0 - x[0], 0.30276E-07)*pow(pip_z[0], -0.66787)*pow(1.0 - pip_z[0], 3.5868)*pow(log(QQ[0]/(0.25*0.25))/log(1.0/(0.25*0.25)), 6.8777/x[0]);
float Ac = (h3*GTMD*pi_thresh*sqrt(pip_PT2[0]/QQ[0])*(2.0 - y[0])*rtz)/term0;
float Acc = (h4*GTMD*pi_thresh*(pip_PT2[0]/QQ[0])*rtz*rtz)/term0;

if(Ac > 0.45) Ac = 0.45; // don't want huge moments
if(Ac < -0.45) Ac = -0.45;
if(Acc > 0.45) Acc = 0.45;
if(Acc < -0.45) Acc = -0.45;

pip_weight = 1.0 + Ac*cos(pip_phih[0]*pi180) + Acc*cos(2.0*pip_phih[0]*pi180);
if(pip_weight != pip_weight) pip_weight = 1.0; // to remove the occasional nan
}
term0file.close();
}

// iteration 2 = from Nick
// pip:
if(ExpOrSim == 0 && accIterationN == 2 && xBin[0] > -1 && QQBin[0] > -1 && W[0] > WMin && y[0] < yMax && pip_index[0] >= 0)
{
ifstream term0file(Form("/scratch/dflt_term0_fromHAPRAD/dflt_term0_fromHAPRAD_BiSc%i_x%iQQ%iz%iPT2%i.txt", binSchemeOpt, xBin[0], QQBin[0], pip_zBin[0], pip_PT2Bin[0]));
if(term0file)
{
float term0;
term0file>>term0;

// from semi_inclusive_model.f: "Parton transverse momentum distribution (Gaussian fit)"
float ac = 1.2025E-10;
float bc = -5.2703E-02;
float cc = 3.7467E-01;
float dc = 6.5397E-02;
float ec = -2.2136E-01;
float fc = -1.0621E-01;

float sgmpt = ac + bc*x[0] + cc*pip_z[0] + dc*x[0]*x[0] + ec*pip_z[0]*pip_z[0] + fc*x[0]*pip_z[0];
if(sgmpt < 0.02) sgmpt = 0.02;
if(sgmpt > 0.15) sgmpt = 0.15;
float GTMD = exp((-(pip_PT2[0] + 0.0))/(2.0*sgmpt))*(1.0/(2.0*pi*sgmpt)); // the 0.0 part comes from "2*pl^2" in haprad, but "pl" seems to always be zero
float pi_thresh = sqrt(1.0 - ((prot_mass + pip_mass)*(prot_mass + pip_mass))/pow(MX_epipX[0], 2.0));
float rt = 1.0 - y[0] - (prot_mass*x[0]*y[0])/(2.0*Beam_Energy);
float rtz = sqrt(rt/(1.0 + ((2.0*prot_mass*x[0])/(y[0]*Beam_Energy))));
float h3 = -0.36544E-03*pow(x[0], -2.1855)*pow(1.0 - x[0], 3.4176)*pow(pip_z[0], -1.7567)*pow(1.0 - pip_z[0], 1.1272)*pow(log(QQ[0]/(0.25*0.25))/log(1.0/(0.25*0.25)), 8.9985);
h3 = h3*12.0*(-0.03 + x[0]); // from Nick: H3_Updated = h3(default)* 12*(-0.03 + x)
float h4 = 0.10908E-02*pow(x[0], -0.35265E-06)*pow(1.0 - x[0], 0.30276E-07)*pow(pip_z[0], -0.66787)*pow(1.0 - pip_z[0], 3.5868)*pow(log(QQ[0]/(0.25*0.25))/log(1.0/(0.25*0.25)), 5.1777/x[0]);
h4 = h4*6.0*(2.0 - x[0]); // from Nick: H4_Updated = 6*(2-x)*h4[default, bb = 5.1777] ... note that bb is changed above
float Ac = (h3*GTMD*pi_thresh*sqrt(pip_PT2[0]/QQ[0])*(2.0 - y[0])*rtz)/term0;
float Acc = (h4*GTMD*pi_thresh*(pip_PT2[0]/QQ[0])*rtz*rtz)/term0;

if(Ac > 0.45) Ac = 0.45; // don't want huge moments
if(Ac < -0.45) Ac = -0.45;
if(Acc > 0.45) Acc = 0.45;
if(Acc < -0.45) Acc = -0.45;

pip_weight = 1.0 + Ac*cos(pip_phih[0]*pi180) + Acc*cos(2.0*pip_phih[0]*pi180);
if(pip_weight != pip_weight) pip_weight = 1.0; // to remove the occasional nan
}
term0file.close();
}

// iteration 1 = default model in haprad:
// pim:
if(ExpOrSim == 0 && accIterationN == 1 && xBin[0] > -1 && QQBin[0] > -1 && W[0] > WMin && y[0] < yMax && pim_index[0] >= 0)
{
ifstream term0file(Form("/scratch/dflt_term0_fromHAPRAD/dflt_term0_fromHAPRAD_BiSc%i_x%iQQ%iz%iPT2%i.txt", binSchemeOpt, xBin[0], QQBin[0], pim_zBin[0], pim_PT2Bin[0]));
if(term0file)
{
float term0;
term0file>>term0;

// from semi_inclusive_model.f: "Parton transverse momentum distribution (Gaussian fit)"
float ac = 1.2025E-10;
float bc = -5.2703E-02;
float cc = 3.7467E-01;
float dc = 6.5397E-02;
float ec = -2.2136E-01;
float fc = -1.0621E-01;

float sgmpt = ac + bc*x[0] + cc*pim_z[0] + dc*x[0]*x[0] + ec*pim_z[0]*pim_z[0] + fc*x[0]*pim_z[0];
if(sgmpt < 0.02) sgmpt = 0.02;
if(sgmpt > 0.15) sgmpt = 0.15;
float GTMD = exp((-(pim_PT2[0] + 0.0))/(2.0*sgmpt))*(1.0/(2.0*pi*sgmpt)); // the 0.0 part comes from "2*pl^2" in haprad, but "pl" seems to always be zero
float pi_thresh = sqrt(1.0 - ((prot_mass + pim_mass)*(prot_mass + pim_mass))/pow(MX_epimX[0], 2.0));
float rt = 1.0 - y[0] - (prot_mass*x[0]*y[0])/(2.0*Beam_Energy);
float rtz = sqrt(rt/(1.0 + ((2.0*prot_mass*x[0])/(y[0]*Beam_Energy))));
float h3 = -0.36544E-03*pow(x[0], -2.1855)*pow(1.0 - x[0], 3.4176)*pow(pim_z[0], -1.7567)*pow(1.0 - pim_z[0], 1.1272)*pow(log(QQ[0]/(0.25*0.25))/log(1.0/(0.25*0.25)), 8.9985);
float h4 = 0.10908E-02*pow(x[0], -0.35265E-06)*pow(1.0 - x[0], 0.30276E-07)*pow(pim_z[0], -0.66787)*pow(1.0 - pim_z[0], 3.5868)*pow(log(QQ[0]/(0.25*0.25))/log(1.0/(0.25*0.25)), 6.8777/x[0]);
float Ac = (h3*GTMD*pi_thresh*sqrt(pim_PT2[0]/QQ[0])*(2.0 - y[0])*rtz)/term0;
float Acc = (h4*GTMD*pi_thresh*(pim_PT2[0]/QQ[0])*rtz*rtz)/term0;

if(Ac > 0.45) Ac = 0.45; // don't want huge moments
if(Ac < -0.45) Ac = -0.45;
if(Acc > 0.45) Acc = 0.45;
if(Acc < -0.45) Acc = -0.45;

pim_weight = 1.0 + Ac*cos(pim_phih[0]*pi180) + Acc*cos(2.0*pim_phih[0]*pi180);
if(pim_weight != pim_weight) pim_weight = 1.0; // to remove the occasional nan
}
term0file.close();
}

// iteration 2 = from Nick
// pim:
if(ExpOrSim == 0 && accIterationN == 2 && xBin[0] > -1 && QQBin[0] > -1 && W[0] > WMin && y[0] < yMax && pim_index[0] >= 0)
{
ifstream term0file(Form("/scratch/dflt_term0_fromHAPRAD/dflt_term0_fromHAPRAD_BiSc%i_x%iQQ%iz%iPT2%i.txt", binSchemeOpt, xBin[0], QQBin[0], pim_zBin[0], pim_PT2Bin[0]));
if(term0file)
{
float term0;
term0file>>term0;

// from semi_inclusive_model.f: "Parton transverse momentum distribution (Gaussian fit)"
float ac = 1.2025E-10;
float bc = -5.2703E-02;
float cc = 3.7467E-01;
float dc = 6.5397E-02;
float ec = -2.2136E-01;
float fc = -1.0621E-01;

float sgmpt = ac + bc*x[0] + cc*pim_z[0] + dc*x[0]*x[0] + ec*pim_z[0]*pim_z[0] + fc*x[0]*pim_z[0];
if(sgmpt < 0.02) sgmpt = 0.02;
if(sgmpt > 0.15) sgmpt = 0.15;
float GTMD = exp((-(pim_PT2[0] + 0.0))/(2.0*sgmpt))*(1.0/(2.0*pi*sgmpt)); // the 0.0 part comes from "2*pl^2" in haprad, but "pl" seems to always be zero
float pi_thresh = sqrt(1.0 - ((prot_mass + pim_mass)*(prot_mass + pim_mass))/pow(MX_epimX[0], 2.0));
float rt = 1.0 - y[0] - (prot_mass*x[0]*y[0])/(2.0*Beam_Energy);
float rtz = sqrt(rt/(1.0 + ((2.0*prot_mass*x[0])/(y[0]*Beam_Energy))));
float h3 = -0.36544E-03*pow(x[0], -2.1855)*pow(1.0 - x[0], 3.4176)*pow(pim_z[0], -1.7567)*pow(1.0 - pim_z[0], 1.1272)*pow(log(QQ[0]/(0.25*0.25))/log(1.0/(0.25*0.25)), 8.9985);
// from Nick: For H3 default values look reasonable.
float h4 = 0.10908E-02*pow(x[0], -0.35265E-06)*pow(1.0 - x[0], 0.30276E-07)*pow(pim_z[0], -0.66787)*pow(1.0 - pim_z[0], 3.5868)*pow(log(QQ[0]/(0.25*0.25))/log(1.0/(0.25*0.25)), 2.077/x[0]);
h4 = h4*75.0*(2.0 - x[0]); // from Nick: H4_updated = 75*(2-x)*h4_default[bb = 2.077]
float Ac = (h3*GTMD*pi_thresh*sqrt(pim_PT2[0]/QQ[0])*(2.0 - y[0])*rtz)/term0;
float Acc = (h4*GTMD*pi_thresh*(pim_PT2[0]/QQ[0])*rtz*rtz)/term0;

if(Ac > 0.45) Ac = 0.45; // don't want huge moments
if(Ac < -0.45) Ac = -0.45;
if(Acc > 0.45) Acc = 0.45;
if(Acc < -0.45) Acc = -0.45;

pim_weight = 1.0 + Ac*cos(pim_phih[0]*pi180) + Acc*cos(2.0*pim_phih[0]*pi180);
if(pim_weight != pim_weight) pim_weight = 1.0; // to remove the occasional nan
}
term0file.close();
}


// phih histos:
for(int j = 1; j >= ExpOrSim; j--)
{
if(xBin[j] > -1 && QQBin[j] > -1 && pip_zBin[j] > -1 && pip_PT2Bin[j] > -1) h_pip_phih[xBin[j]][QQBin[j]][pip_zBin[j]][pip_PT2Bin[j]][j]->Fill(pip_phih[j], pip_weight);
if(xBin[j] > -1 && QQBin[j] > -1 && pim_zBin[j] > -1 && pim_PT2Bin[j] > -1) h_pim_phih[xBin[j]][QQBin[j]][pim_zBin[j]][pim_PT2Bin[j]][j]->Fill(pim_phih[j], pim_weight);
}

// %%%%%%%%% end analysis of kinematics %%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(i%1000 == 0) cout<<"\ranalyzed "<<i/1000<<" thousand events"<<flush;
} // end if(genExclusiveEvent == 0)
} // end of loop over entries
// %%%%%%%%%%%%%%%%%%%%%%%%
// %%%%% END POTATOES %%%%%

outputfile->Write();
cout<<endl<<"Done!"<<endl;
stopwat->Print();

} // end of program
