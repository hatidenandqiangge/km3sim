#ifndef KM3EMEnergyFlux_h
#define KM3EMEnergyFlux_h 1

#include <vector>
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "KM3EMDistanceFlux.hh"

class KM3EMEnergyFlux {
public:
  KM3EMEnergyFlux(char *, double, double, int, double);
  ~KM3EMEnergyFlux();

public:
  void FindBins(double energyin, double distancein, double anglein);
  int GetNumberOfSamples() { return NumberOfSamples; };
  onePE GetSamplePoint();
  bool ModelTrigger(double TheE);

private:
  int NEnergies;
  std::vector<KM3EMDistanceFlux *> *keepEnergies;
  int ibin1;
  int ibin2;
  double Flux;
  double FluxRMS;
  double ratio;
  double EnergyMin;
  double EnergyMax;
  int NumberOfSamples;
  double RatioThis;
  double MaxBoostRatio;
};

#endif
