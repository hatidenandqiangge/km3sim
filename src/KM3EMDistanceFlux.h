#ifndef KM3EMDistanceFlux_h
#define KM3EMDistanceFlux_h 1

#include <vector>
#include "globals.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "KM3EMAngularFlux.h"

class KM3EMDistanceFlux {
 public:
  KM3EMDistanceFlux(std::ifstream &);
  ~KM3EMDistanceFlux();

 public:
  void FindBins(G4double distancein, G4double anglein);
  G4double GiveEnergy() { return Energy; };
  G4double GiveFlux() { return Flux; };
  G4double GiveFluxRMS() { return FluxRMS; };
  onePE GetSamplePoint();

 private:
  std::vector<KM3EMAngularFlux *> *keepDistances;
  G4double Energy;
  G4int ibin1;
  G4int ibin2;
  G4double Flux;
  G4double FluxRMS;
  G4double ratio;
  G4int VertexDistanceBins;
};

#endif
