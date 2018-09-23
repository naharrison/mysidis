
// returns the index of the first reconstructed particle that closely matches the provided generated particle
Int_t getIndexOfMatch(Int_t gq, Float_t gp, Float_t gth, Float_t gph, Int_t gpart, Int_t q[], Float_t p[], Float_t cx[], Float_t cy[], Float_t cz[]) {
  TVector3 g3 = TVector3(1.0, 1.0, 1.0);
  g3.SetTheta((3.14159/180.0)*gth);
  g3.SetPhi((3.14159/180.0)*gph);
  g3.SetMag(gp);
  int r_index = -123;
  for(int ir = 0; ir < gpart; ir++) {
    if(q[ir] == gq) {
      TVector3 r3 = TVector3(p[ir]*cx[ir], p[ir]*cy[ir],p[ir]*cz[ir]); // candidate
      if(fabs(g3.Mag() - r3.Mag()) < 1.5*0.04 && fabs(g3.Theta()-r3.Theta()) < 1.5*0.007 && fabs(g3.Phi()-r3.Phi()) < 1.5*0.05) {
        r_index = ir;
        break;
      }
    }
  }

  return r_index;
}



Int_t getIndexOfPid(Int_t pid, Int_t mcnentr, Int_t mcid[]) {
  int index = -123;
  for(int ig = 0; ig < mcnentr; ig++) {
    if(mcid[ig] == pid) {
      index = ig;
      break;
    }
  }

  return index;
}



void printInfo(Int_t eIndex, Int_t hadIndex, Int_t hadID, Int_t gpart, Float_t p[], Float_t cx[], Float_t cy[], Float_t cz[], Float_t sc_t[], Float_t sc_r[], Float_t ec_ei[], Float_t ec_eo[], UShort_t nphe[]) {

  float speed_of_light = 29.9792458; // cm/ns
  float startTime = sc_t[eIndex] - (sc_r[eIndex]/speed_of_light);
  float hadTOF = sc_t[hadIndex] - startTime;
  float velocity = (sc_r[hadIndex]/hadTOF)/speed_of_light;

  TVector3 had3 = TVector3(p[hadIndex]*cx[hadIndex], p[hadIndex]*cy[hadIndex],p[hadIndex]*cz[hadIndex]);

  // if statements are loose cuts to remove a small number of unrealistic outliers
  if(p[hadIndex] > 0.21 && p[hadIndex] < 5.3 && had3.Theta() > 0.088 && had3.Theta() < 2.22 && velocity > 0 && velocity < 1.5) {
    if(nphe[hadIndex] >= 0 && nphe[hadIndex] < 350 && ec_ei[hadIndex] >= 0 && ec_ei[hadIndex] < 0.9 && ec_eo[hadIndex] >= 0 && ec_eo[hadIndex] < 1.1) {
      cout << hadID << "," << p[hadIndex] << "," << had3.Theta() << "," << velocity << "," << nphe[hadIndex] << "," << ec_ei[hadIndex] << "," << ec_eo[hadIndex] << endl;
    }
  }

}
