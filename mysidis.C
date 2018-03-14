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
int Nfiles = 5;
int ExpOrSim = 0; // 0(MC) or 1(data)

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

bool genExclusiveEvent = 0; // if the generated event is an exclusive event, i.e. ep --> e pi+ n, we don't want to analyze it at all (the generator give unrealistic kinematics)
if(ExpOrSim == 0 && mcnentr == 3)
{
if((mcid[1] == 211 && mcid[2] == 2112) || (mcid[1] == 2112 && mcid[2] == 211)) genExclusiveEvent = 1; // the first (zeroth) generated particle is always the scattered electron
}

if(genExclusiveEvent == 0)
{
cout<<"event #"<<i<<endl;

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

cout<<"e- indices (gen rec): "<<e_index[0]<<" "<<e_index[1]<<endl;

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
}
if(pim_index[1] > -122)
{
V3_pim[1].SetXYZ(p[pim_index[1]]*cx[pim_index[1]], p[pim_index[1]]*cy[pim_index[1]], p[pim_index[1]]*cz[pim_index[1]]);
V4_pim[1].SetXYZT(V3_pim[1].X(), V3_pim[1].Y(), V3_pim[1].Z(), sqrt(V3_pim[1].Mag2() + pow(pim_mass,2)));
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

} // end if(genExclusiveEvent == 0)
cout<<endl;
} // end of loop over entries
} // end of program
