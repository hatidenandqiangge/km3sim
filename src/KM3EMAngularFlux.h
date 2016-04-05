#ifndef KM3EMAngularFlux_h
#define KM3EMAngularFlux_h 1

#include <vector>
#include "globals.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "KM3EMTimePointDis.h"

class KM3EMAngularFlux {
 public:
  KM3EMAngularFlux(std::ifstream &, bool &ok, bool FineBin);
  ~KM3EMAngularFlux();

 public:
  void FindBins(G4double anglein);
  G4double GiveDistance() { return Distance; };
  G4double GiveFlux() { return Flux; };
  G4double GiveFluxRMS() { return FluxRMS; };
  bool IsValid() { return IsThisValid; };
  onePE GetSamplePoint();

 private:
  std::vector<KM3EMTimePointDis *> *keepAngles;
  G4int VertexSolidAngleBins;
  G4double Distance;
  bool IsThisValid;
  G4int ibin1;
  G4int ibin2;
  G4double Flux;
  G4double FluxRMS;
  G4double ratio;
};

#endif
