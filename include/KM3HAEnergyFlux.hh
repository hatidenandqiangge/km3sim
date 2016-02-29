#ifndef KM3HAEnergyFlux_h
#define KM3HAEnergyFlux_h 1

#include <vector>
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "KM3EMDistanceFlux.hh"

class KM3HAEnergyFlux
{
public:
  KM3HAEnergyFlux(char *,G4double,G4double,G4double,G4double);
  ~KM3HAEnergyFlux();
public:
  void FindBins(G4int idbeam, G4double energyin,G4double distancein,G4double anglein);
  G4int GetNumberOfSamples(){return NumberOfSamples;};
  onePE GetSamplePoint();
private:
  G4double Rescale(G4int idbeam, G4double energyin,G4double energydis);
  G4int NPartsDists;
  std::vector<KM3EMDistanceFlux*>* keepEnergies;
  G4double ParticleEnergyScaleFactors[7][5];
  G4int ibin;
  G4double Flux;
  G4double FluxRMS;
  G4double EnergyMin;
  G4double EnergyMax;
  G4int NumberOfSamples;
  G4double RatioThis;
};

#endif
