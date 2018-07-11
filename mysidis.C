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

void mysidis()
{
int Nfiles = 1;

// %%%%% read the input files into a TChain %%%%%
string firstfilename = "";
string lastfilename = "";

ifstream filelist("programFiles/v12_MCFiles.txt");

TChain *h22chain = new TChain("h22");

for(int k = 0; k < Nfiles; k++)
{
string filename;
filelist>>filename;
h22chain->Add(filename.c_str());
}
// %%%%% end read the input files into a TChain %%%%%

// %%%%% define variables %%%%%
Float_t e_mass = 0.000511; // GeV
Float_t prot_mass = 0.938272; // GeV
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
// %%%%% end define variables %%%%%

// %%%%% set branch addresses %%%%%
h22chain->SetBranchAddress("mcnentr", &mcnentr);
h22chain->SetBranchAddress("mcid", mcid);
h22chain->SetBranchAddress("mctheta", mctheta);
h22chain->SetBranchAddress("mcphi", mcphi);
h22chain->SetBranchAddress("mcp", mcp);
// %%%%% end set branch addresses %%%%%

// %%%%% define histograms %%%%%
TH2F *hQQvsx = new TH2F("hQQvsx", "hQQvsx", 100, 0, 1, 100, 0, 5);
// %%%%% end define histograms %%%%%

for(int i = 0; i < h22chain->GetEntries(); i++) // loop over entries
{
h22chain->GetEntry(i);

bool genExclusiveEvent = 0; // if the generated event is an exclusive event, i.e. ep --> e pi+ n, we don't want to analyze it at all (the generator give unrealistic kinematics)
if(mcnentr == 3)
{
if((mcid[1] == 211 && mcid[2] == 2112) || (mcid[1] == 2112 && mcid[2] == 211)) genExclusiveEvent = 1; // the first (zeroth) generated particle is always the scattered electron
}

if(genExclusiveEvent == 0)
{

TLorentzVector V4k(0.0, 0.0, Beam_Energy, Beam_Energy); // x, y, z, t
TLorentzVector V4ISproton(0.0, 0.0, 0.0, prot_mass); // IS = Initial State

// %%%%% set e- 3 and 4 vectors %%%%%
TVector3 V3_e;
TLorentzVector V4_e;

V3_e.SetXYZ(mcp[0]*cos(pi180*mcphi[0])*sin(pi180*mctheta[0]), mcp[0]*sin(pi180*mcphi[0])*sin(pi180*mctheta[0]), mcp[0]*cos(pi180*mctheta[0]));
V4_e.SetXYZT(V3_e.X(), V3_e.Y(), V3_e.Z(), sqrt(V3_e.Mag2() + pow(e_mass, 2)));

TVector3 V3_H;
TLorentzVector V4_q, V4_H;

V4_q = V4k - V4_e;
V4_H = V4_q + V4ISproton;
V3_H.SetXYZ(V4_H.X(), V4_H.Y(), V4_H.Z());
// %%%%% end set e- 3 and 4 vectors %%%%%

// %%%%%%%%% analysis of kinematics %%%%%%%%%%%
float W = V4_H.Mag();
float y = V4_q.T()/Beam_Energy;
float QQ = -V4_q*V4_q;
float x = QQ/(2.0*V4_q.T()*prot_mass);
hQQvsx->Fill(x, QQ);
// %%%%%%%%% end analysis of kinematics %%%%%%%%%

} // end if(genExclusiveEvent == 0)
} // end of loop over entries

hQQvsx->Draw("colz");

} // end of program
