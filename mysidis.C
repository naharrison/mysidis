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
#include "programFiles/eID.C"
#include "programFiles/hadronID.C"
#include "programFiles/getGenIndices.C"

void mysidis(int e_zvertex_strict = 0, int e_ECsampling_strict = 0, int e_ECoutVin_strict = 0, int e_ECgeometric_strict = 0, int e_CCthetaMatching_strict = 0, int e_R1fid_strict = 0, int e_R3fid_strict = 0, int e_CCphiMatching_strict = 0, int e_CCfiducial_strict = 0, int yCut_strict = 0, int pip_vvp_strict = 0, int pip_R1fid_strict = 0, int pip_MXcut_strict = 0, int pim_vvp_strict = 0, int pim_R1fid_strict = 0, int pim_MXcut_strict = 0)
{
int MCversion = 12;
int filestart = 1;
int Nfiles = 1;
int ExpOrSim = 0; // 0(MC) or 1(data)

// %%%%% define some historams %%%%%
TH1F *hpres = new TH1F("hpres", "hpres", 100, -0.1, 0.1);
TH1F *hthetares = new TH1F("hthetares", "hthetares", 100, -0.02, 0.02);
TH1F *hphires = new TH1F("hphires", "hphires", 100, -0.1, 0.1);
// %%%% end define historams %%%%%%%

// %%%%% some cut values %%%%%
float WMin = 2.05;
float yMax = 0.85; // default value (for when yCut_strict = 0)
if(yCut_strict == 1) yMax = 0.8;
// %%%%% end some cut values %%%%%

// %%%%% read the input files into a TChain %%%%%
int NtotalFiles = 5;

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

  cout << "event #" << i << endl;

  cout << "  generated particles:" << endl;
  for(int ig = 0; ig < mcnentr; ig++) {
    TVector3 vec3 = TVector3(1.0, 1.0, 1.0);
    vec3.SetTheta((3.14159/180.0)*mctheta[ig]);
    vec3.SetPhi((3.14159/180.0)*mcphi[ig]);
    vec3.SetMag(mcp[ig]);
    cout << "    " << mcid[ig] << ": p=" << vec3.Mag() << ", theta=" << vec3.Theta() << ", phi=" << vec3.Phi() << endl;
  }


  cout << endl << "  reconstructed tracks:" << endl;
  for(int ir = 0; ir < gpart; ir++) {
    TVector3 vec3 = TVector3(p[ir]*cx[ir], p[ir]*cy[ir],p[ir]*cz[ir]);
    cout << "    q=" << q[ir] << ", p=" << vec3.Mag() << ", theta=" << vec3.Theta() << ", phi=" << vec3.Phi() << endl;
  }


  cout << endl << "  comparison of gen and rec:" << endl;
  for(int ig = 0; ig < mcnentr; ig++) {
    for(int ir = 0; ir < gpart; ir++) {
      cout << "    g" << ig << " r" << ir << " comparison: ";

      TVector3 rvec3 = TVector3(p[ir]*cx[ir], p[ir]*cy[ir],p[ir]*cz[ir]);
      TVector3 gvec3 = TVector3(1.0, 1.0, 1.0);
      gvec3.SetTheta((3.14159/180.0)*mctheta[ig]);
      gvec3.SetPhi((3.14159/180.0)*mcphi[ig]);
      gvec3.SetMag(mcp[ig]);

      if(fabs(gvec3.Mag() - rvec3.Mag()) < 0.04 && fabs(gvec3.Theta()-rvec3.Theta()) < 0.007 && fabs(gvec3.Phi()-rvec3.Phi()) < 0.05) cout << "the momenta are close" << endl;
      else cout << "the momenta are not close" << endl;

      hpres->Fill(gvec3.Mag() - rvec3.Mag());
      hthetares->Fill(gvec3.Theta() - rvec3.Theta());
      hphires->Fill(gvec3.Phi() - rvec3.Phi());
    }
  }

  
  cout << endl;

} // end of loop over entries


hpres->Draw();
new TCanvas();
hthetares->Draw();
new TCanvas();
hphires->Draw();


} // end of program
