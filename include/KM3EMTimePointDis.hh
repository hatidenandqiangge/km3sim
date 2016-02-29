#ifndef KM3EMTimePointDis_h
#define KM3EMTimePointDis_h 1

#include <vector>
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "G4ExceptionHandler.hh"

struct onePE{
  G4double time;
  G4double costh;
  G4double phi;
};

class KM3EMTimePointDis
{
public:
  KM3EMTimePointDis(std::ifstream&, bool& ok);
  ~KM3EMTimePointDis();
public:
  onePE GetSamplePoint();
  G4double GiveAngle(){return angle;};
  G4double GiveFlux(){return Flux;};
  G4double GiveFluxRMS(){return FluxRMS;};
  bool IsValid(){return IsThisValid;};
private:
  std::vector<G4float>* keepDis;
  std::vector<G4float>* keepTh2Th3Num;
  std::vector<G4float>* keepExpoTh2;
  std::vector<G4float>* keepExpoTh3;
  bool time_ok[52];
  G4double pi2;
  G4double angle;
  G4double Flux;
  G4double FluxRMS;
  bool IsThisValid;
  //definition of bins limits for direction sampling
  G4double theta_Low[834]; //this should be OMSolidAngleBins below
  G4double theta_High[834];
  G4double phi_Low[834];
  G4double phi_High[834];
  G4int TimeSolidAngleBins;
  G4int TimeBins;
  G4int OMSolidAngleBins;
  G4int TimeTimeSolidAngleBins;
  //////////////////////////////////////////////////
};

#endif
