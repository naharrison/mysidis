#include <fstream>
#include <iostream>
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
#include "programFiles/pidFunctions.C"

void mysidis()
{
int MCversion = 12;
int Nfiles = 100;
int ExpOrSim = 0; // 0(MC) or 1(data)

// %%%%% read the input files into a TChain %%%%%
ifstream filelist;
if(ExpOrSim == 0) filelist.open(Form("programFiles/v%i_MCFiles.txt", MCversion));
else filelist.open("programFiles/dataFiles.txt");

TChain *h22chain = new TChain("h22");

for(int k = 0; k < Nfiles; k++) {
  string filename;
  filelist>>filename;
  h22chain->Add(filename.c_str());
}
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

Int_t mcnentr, mcid[35];
Float_t mctheta[35], mcphi[35], mcp[35];

Int_t gpart;

Float_t p[35];
Int_t q[35];
Float_t cx[35],  cy[35],  cz[35],  vz[35],  vy[35],  vx[35],  b[35];

UChar_t dc_sect[35];
Float_t tl1_x[35], tl1_y[35];
Float_t tl3_x[35]; // used to be dc_xsc in h10
Float_t tl3_y[35], tl3_z[35];
Float_t tl3_cx[35]; // used to be dc_cxsc in h10
Float_t tl3_cy[35], tl3_cz[35];

UChar_t ec_sect[35];
Float_t etot[35],  ec_ei[35],  ec_eo[35],  ec_t[35],  ech_x[35],  ech_y[35],  ech_z[35];

UChar_t sc_sect[35];
Float_t sc_t[35], sc_r[35];
UChar_t sc_pd[35];

UChar_t cc_sect[35];
UShort_t cc_segm[35], nphe[35];
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

// %%%%% POTATOES %%%%%
// %%%%%%%%%%%%%%%%%%%%
for(int i = 0; i < h22chain->GetEntries(); i++) // loop over entries
{
  h22chain->GetEntry(i);

  int genPip_i = getIndexOfPid(211, mcnentr, mcid);
  int genKp_i = getIndexOfPid(321, mcnentr, mcid);
  int genP_i = getIndexOfPid(2212, mcnentr, mcid);
  int genPos_i = getIndexOfPid(-11, mcnentr, mcid);

  int recE_i = getIndexOfMatch(-1, mcp[0], mctheta[0], mcphi[0], gpart, q, p, cx, cy, cz);

  int recPip_i = -123;
  if(genPip_i > -1) recPip_i = getIndexOfMatch(1, mcp[genPip_i], mctheta[genPip_i], mcphi[genPip_i], gpart, q, p, cx, cy, cz);

  int recKp_i = -123;
  if(genKp_i > -1) recKp_i = getIndexOfMatch(1, mcp[genKp_i], mctheta[genKp_i], mcphi[genKp_i], gpart, q, p, cx, cy, cz);

  int recP_i = -123;
  if(genP_i > -1) recP_i = getIndexOfMatch(1, mcp[genP_i], mctheta[genP_i], mcphi[genP_i], gpart, q, p, cx, cy, cz);

  int recPos_i = -123;
  if(genPos_i > -1) recPos_i = getIndexOfMatch(1, mcp[genPos_i], mctheta[genPos_i], mcphi[genPos_i], gpart, q, p, cx, cy, cz);


  if(recE_i > -1 && recPip_i > -1) {
    printInfo(recE_i, recPip_i, 211, gpart, p, cx, cy, cz, sc_t, sc_r, ec_ei, ec_eo, nphe);
  }

  if(recE_i > -1 && recKp_i > -1) {
    printInfo(recE_i, recKp_i, 321, gpart, p, cx, cy, cz, sc_t, sc_r, ec_ei, ec_eo, nphe);
  }

  if(recE_i > -1 && recP_i > -1) {
    printInfo(recE_i, recP_i, 2212, gpart, p, cx, cy, cz, sc_t, sc_r, ec_ei, ec_eo, nphe);
  }

  if(recE_i > -1 && recPos_i > -1) {
    printInfo(recE_i, recPos_i, -11, gpart, p, cx, cy, cz, sc_t, sc_r, ec_ei, ec_eo, nphe);
  }



} // end of loop over entries


} // end of program
