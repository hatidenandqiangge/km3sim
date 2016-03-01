#ifndef KM3EMDirectFlux_h
#define KM3EMDirectFlux_h 1

#include <vector>
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "KM3EMAngularFlux.hh"

class KM3EMDirectFlux {
public:
  KM3EMDirectFlux(char *, G4double);
  ~KM3EMDirectFlux();

public:
  void FindBins(G4double MeanNumPhotons, G4double distancein, G4double anglein);
  G4int GetNumberOfSamples() { return NumberOfSamples; };
  onePE GetSamplePoint();

private:
  std::vector<KM3EMAngularFlux *> *keepDistances;
  G4int ibin1;
  G4int ibin2;
  G4double Flux;
  G4double ratio;
  G4int NumberOfSamples;
  G4double RatioThis;
  G4int VertexDistanceBins;
};

#endif
