#include "KM3EMTimePointDis.h"
#include "Randomize.hh"

using CLHEP::degree;
using CLHEP::pi;

KM3EMTimePointDis::KM3EMTimePointDis(std::ifstream &infile, bool &ok) {
  TimeSolidAngleBins = 52;
  TimeBins = 111;
  OMSolidAngleBins = 834;
  TimeTimeSolidAngleBins = TimeSolidAngleBins * TimeBins;
  keepDis = new std::vector<G4float>;
  keepDis->reserve(TimeTimeSolidAngleBins);
  keepTh2Th3Num = new std::vector<G4float>;
  keepTh2Th3Num->reserve(OMSolidAngleBins);
  keepExpoTh2 = new std::vector<G4float>;
  keepExpoTh2->reserve(OMSolidAngleBins);
  keepExpoTh3 = new std::vector<G4float>;
  keepExpoTh3->reserve(OMSolidAngleBins);

  G4float val;
  char valC[4];
  infile.read(valC, 4);
  angle = double(*(float *)valC);
  infile.read(valC, 4);
  Flux = double(*(float *)valC);
  infile.read(valC, 4);
  FluxRMS = double(*(float *)valC);
  if (FluxRMS <= 0.0) FluxRMS = sqrt(Flux);
  FluxRMS = FluxRMS * FluxRMS;
  Flux = log(Flux);
  FluxRMS = log(FluxRMS);
  for (G4int i = 0; i < OMSolidAngleBins; i++) {
    infile.read(valC, 4);
    val = *(float *)valC;
    keepTh2Th3Num->push_back(val);
    infile.read(valC, 4);
    val = *(float *)valC;
    keepExpoTh2->push_back(val);
    infile.read(valC, 4);
    val = *(float *)valC;
    keepExpoTh3->push_back(val);
  }
  for (G4int i = 0; i < TimeTimeSolidAngleBins; i++) {
    infile.read(valC, 4);
    val = *(float *)valC;
    keepDis->push_back(val);
  }
  pi2 = 2.0 * pi;

  ok = ((*keepTh2Th3Num)[OMSolidAngleBins - 1] > 0.999);
  if (ok)
    IsThisValid = true;
  else {
    IsThisValid = false;
    keepTh2Th3Num->clear();
    delete keepTh2Th3Num;
    keepTh2Th3Num = NULL;
    keepDis->clear();
    delete keepDis;
    keepDis = NULL;
  }

  for (G4int it23 = 0; it23 < TimeSolidAngleBins; it23++) {
    G4int ilast = it23 * TimeBins + TimeBins - 1;
    time_ok[it23] = ((*keepDis)[ilast] > 0.999);
    if (!time_ok[it23]) G4cout << "Time bin is null " << it23 << G4endl;
  }
  //  G4cout<<"KM3EMTimePointDis "<<angle<<" "<<Flux<<"
  //  "<<(*keepTh2Th3Num)[833]<<G4endl;
  // definition of bins limits for direction sampling
  // definition of two solid angle areas, one with 3 degrees binning and the
  // second with 6 degrees binning
  G4int NumberOfThetas1 = 34;
  G4int NumberOfThetas2 = 13;
  G4double Theta1Min = 0.0;
  G4double Theta1Max = 102.0;
  G4double Theta2Min = 102.0;
  G4double Theta2Max = 180.0;

  G4int NumberOfThetas = NumberOfThetas1 + NumberOfThetas2;
  G4int ibinNum_Tot = 0;
  // first solid angle area
  G4double dtheta = (Theta1Max - Theta1Min) / NumberOfThetas1;
  G4double cosdtheta = cos(dtheta * degree);
  for (G4int ith = 0; ith < NumberOfThetas1; ith++) {
    G4double thetalow = Theta1Min + dtheta * ith;
    G4double thetahigh = thetalow + dtheta;
    G4double costhetalow = fabs(cos(thetalow * degree));
    G4double costhetahigh = fabs(cos(thetahigh * degree));
    G4double cosmin;
    if (costhetalow < costhetahigh)
      cosmin = costhetalow;
    else
      cosmin = costhetahigh;
    cosmin = cosmin * cosmin;
    G4double cosdphi = (cosdtheta - cosmin) / (1 - cosmin);
    G4double dphi = acos(cosdphi) / degree;
    if (dphi < 9.0) dphi = 9.0;
    G4int NumberOfPhis = int(ceil(180.0 / dphi));
    dphi = 180.0 / NumberOfPhis;
    for (G4int iph = 0; iph < NumberOfPhis; iph++) {
      G4double philow = dphi * iph;
      G4double phihigh = philow + dphi;
      theta_Low[ibinNum_Tot] = thetalow;
      theta_High[ibinNum_Tot] = thetahigh;
      if (ith == NumberOfThetas1 - 1) theta_High[ibinNum_Tot] = Theta1Max;
      phi_Low[ibinNum_Tot] = philow;
      phi_High[ibinNum_Tot] = phihigh;
      if (iph == NumberOfPhis - 1) phi_High[ibinNum_Tot] = 180.0;
      ibinNum_Tot++;
    }
  }
  // second solid angle area
  dtheta = (Theta2Max - Theta2Min) / NumberOfThetas2;
  cosdtheta = cos(dtheta * degree);
  for (G4int ith = NumberOfThetas1; ith < NumberOfThetas; ith++) {
    G4double thetalow = Theta2Min + dtheta * (ith - NumberOfThetas1);
    G4double thetahigh = thetalow + dtheta;
    G4double costhetalow = fabs(cos(thetalow * degree));
    G4double costhetahigh = fabs(cos(thetahigh * degree));
    G4double cosmin;
    if (costhetalow < costhetahigh)
      cosmin = costhetalow;
    else
      cosmin = costhetahigh;
    cosmin = cosmin * cosmin;
    G4double cosdphi = (cosdtheta - cosmin) / (1 - cosmin);
    G4double dphi = acos(cosdphi) / degree;
    if (dphi < 9.0) dphi = 9.0;
    G4int NumberOfPhis = int(ceil(180.0 / dphi));
    dphi = 180.0 / NumberOfPhis;
    for (G4int iph = 0; iph < NumberOfPhis; iph++) {
      G4double philow = dphi * iph;
      G4double phihigh = philow + dphi;
      theta_Low[ibinNum_Tot] = thetalow;
      theta_High[ibinNum_Tot] = thetahigh;
      if (ith == NumberOfThetas - 1) theta_High[ibinNum_Tot] = Theta2Max;
      phi_Low[ibinNum_Tot] = philow;
      phi_High[ibinNum_Tot] = phihigh;
      if (iph == NumberOfPhis - 1) phi_High[ibinNum_Tot] = 180.0;
      ibinNum_Tot++;
    }
  }
  if (ibinNum_Tot != OMSolidAngleBins)
    G4Exception(
        "Error calculated direction bins are not the same as in KM3SD\n", "",
        FatalException, "");  // number needs to change
  //////////////////////////////////////////////////
}
KM3EMTimePointDis::~KM3EMTimePointDis() {
  if (keepDis != NULL) {
    keepDis->clear();
    delete keepDis;
    keepDis = NULL;
  }
  if (keepTh2Th3Num != NULL) {
    keepTh2Th3Num->clear();
    delete keepTh2Th3Num;
    keepTh2Th3Num = NULL;
  }
}

// gives the random values. Sampling based on sorting the pdf
onePE KM3EMTimePointDis::GetSamplePoint() {
  onePE aPE;
  G4double time;
  G4double costh;
  G4double theta, phi;
  if (!IsThisValid)
    G4Exception("Error sampling point and time for null distribution\n", "",
                FatalException, "");
  // first we sample a direction point using the cumulative keepTh2Th3Num
  // [0-833]
  G4double rrr = G4UniformRand();
  G4int ibinNum23;
  for (ibinNum23 = 0; ibinNum23 < OMSolidAngleBins; ibinNum23++) {
    if (rrr < (*keepTh2Th3Num)[ibinNum23]) break;
  }
  // next we sample an exponential theta and phi that belongs to this slice
  // the definitions of the bin limits are calculated in the constructor
  // the same definitions are in KM3SD;
  // first sample a theta
  double minval = theta_Low[ibinNum23];
  double maxval = theta_High[ibinNum23];
  double param = (*keepExpoTh2)[ibinNum23];
  if (fabs(param) > 1.0e-6) {
    double dm = param * (maxval - minval);
    if (dm < 700.0)
      theta =
          minval + (1.0 / param) * log(1.0 + G4UniformRand() * (exp(dm) - 1.0));
    else {
      double rrr = G4UniformRand();
      if (rrr > 0.0)
        theta = maxval + log(rrr) / param;
      else
        theta = minval;
    }
  } else
    theta = minval + G4UniformRand() * (maxval - minval);
  // next sample a phi
  minval = phi_Low[ibinNum23];
  maxval = phi_High[ibinNum23];
  param = (*keepExpoTh3)[ibinNum23];
  if (fabs(param) > 1.0e-6) {
    double dm = param * (maxval - minval);
    if (dm < 700.0)
      phi =
          minval + (1.0 / param) * log(1.0 + G4UniformRand() * (exp(dm) - 1.0));
    else {
      double rrr = G4UniformRand();
      if (rrr > 0.0)
        phi = maxval + log(rrr) / param;
      else
        phi = minval;
    }
  } else
    phi = minval + G4UniformRand() * (maxval - minval);
  //  G4cout<<"ibinNum23= "<<ibinNum23<<" theta= "<<theta<<" phi=
  //  "<<phi<<G4endl;
  // next we must find the time bins that belongs this photon
  // follows the limits of the bins
  static double maxth2[17] = {
      4.43924,  // the definitions of the array and the next one is in KM3SD.cc
      6.27957, 8.10961, 9.93636, 11.4783, 14.0699, 18.1949, 22.3316, 25.8419,
      36.8699, 45.5730, 53.1301, 63.2563, 90.0000, 113.578, 143.130, 180.0};
  static double maxth3[52] = {
      55.0,  105.0, 180.0, 45.0,  110.0, 180.0, 45.0,  110.0, 180.0,
      45.0,  95.0,  180.0, 45.0,  120.0, 180.0, 45.0,  110.0, 180.0,
      40.0,  75.0,  180.0, 30.0,  45.0,  95.0,  180.0, 40.0,  75.0,
      110.0, 180.0, 30.0,  45.0,  85.0,  145.0, 180.0, 40.0,  75.0,
      115.0, 180.0, 40.0,  95.0,  180.0, 45.0,  85.0,  180.0, 40.0,
      110.0, 180.0, 65.0,  120.0, 180.0, 180.0, 180.0};  // dimension 52 should
                                                         // be equal to
                                                         // TimeSolidAngleBins
  static int indexes[18] = {0,  3,  6,  9,  12, 15, 18, 21, 25,
                            29, 34, 38, 41, 44, 47, 50, 51, 52};
  G4int ith;
  for (ith = 0; ith < 17; ith++)
    if (theta < maxth2[ith]) break;
  if (ith > 16) ith = 16;
  G4int iph;
  for (iph = indexes[ith]; iph < indexes[ith + 1]; iph++)
    if (phi < maxth3[iph]) break;
  if (iph > indexes[ith + 1] - 1) iph = indexes[ith + 1] - 1;
  //  G4cout<<" iph= "<<iph<<G4endl;
  ////here we must check that the iph time bin in not null and find another
  ////////////////////////////

  G4int cang23bin = iph;
  rrr = G4UniformRand();
  G4int ibint;
  for (ibint = TimeBins * cang23bin; ibint < TimeBins * cang23bin + TimeBins;
       ibint++) {
    if (rrr < (*keepDis)[ibint]) break;
  }
  ibint -= cang23bin * TimeBins;
  if (ibint < 40)
    time = double(ibint) * 0.5 - 10.0 + G4UniformRand() * 0.5;
  else if (ibint < 60)
    time = double(ibint - 40) + 10.0 + G4UniformRand();
  else if (ibint < 80)
    time = double(ibint - 60) * 3.0 + 30.0 + G4UniformRand() * 3.0;
  else if (ibint < 90)
    time = double(ibint - 80) * 11.0 + 90.0 + G4UniformRand() * 11.0;
  else if (ibint < 106)
    time = double(ibint - 90) * 50.0 + 200.0 + G4UniformRand() * 50.0;
  else
    time = double(ibint - 106) * 200.0 + 1000.0 + G4UniformRand() * 200.0;

  // up to here the th2 and th3 is in degrees
  costh = cos(theta * degree);
  phi *= degree;

  if (G4UniformRand() < 0.5) phi = pi2 - phi;
  aPE.time = time;
  aPE.costh = costh;
  aPE.phi = phi;
  return aPE;
}
