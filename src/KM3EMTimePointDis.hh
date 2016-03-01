#ifndef KM3EMTimePointDis_h
#define KM3EMTimePointDis_h 1

#include <vector>
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "G4ExceptionHandler.hh"

struct onePE {
  double time;
  double costh;
  double phi;
};

class KM3EMTimePointDis {
public:
  KM3EMTimePointDis(std::ifstream &, bool &ok);
  ~KM3EMTimePointDis();

public:
  onePE GetSamplePoint();
  double GiveAngle() { return angle; };
  double GiveFlux() { return Flux; };
  double GiveFluxRMS() { return FluxRMS; };
  bool IsValid() { return IsThisValid; };

private:
  std::vector<float> *keepDis;
  std::vector<float> *keepTh2Th3Num;
  std::vector<float> *keepExpoTh2;
  std::vector<float> *keepExpoTh3;
  bool time_ok[52];
  double pi2;
  double angle;
  double Flux;
  double FluxRMS;
  bool IsThisValid;
  // definition of bins limits for direction sampling
  double theta_Low[834]; // this should be OMSolidAngleBins below
  double theta_High[834];
  double phi_Low[834];
  double phi_High[834];
  int TimeSolidAngleBins;
  int TimeBins;
  int OMSolidAngleBins;
  int TimeTimeSolidAngleBins;
  //////////////////////////////////////////////////
};

#endif
