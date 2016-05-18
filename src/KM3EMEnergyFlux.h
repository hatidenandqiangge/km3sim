#ifndef KM3EMEnergyFlux_h
#define KM3EMEnergyFlux_h 1

#include <vector>
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "KM3EMDistanceFlux.h"

class KM3EMEnergyFlux {
 public:
  KM3EMEnergyFlux(char *, G4double, G4double, G4int, G4double);
  ~KM3EMEnergyFlux();

 public:
  void FindBins(G4double energyin, G4double distancein, G4double anglein);
  G4int GetNumberOfSamples() { return NumberOfSamples; };
  onePE GetSamplePoint();
  G4bool ModelTrigger(G4double TheE);

 private:
  G4int NEnergies;
  std::vector<KM3EMDistanceFlux *> *keepEnergies;
  G4int ibin1;
  G4int ibin2;
  G4double Flux;
  G4double FluxRMS;
  G4double ratio;
  G4double EnergyMin;
  G4double EnergyMax;
  G4int NumberOfSamples;
  G4double RatioThis;
  G4double MaxBoostRatio;
};

#endif
