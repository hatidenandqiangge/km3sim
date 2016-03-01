#ifndef KM3HAEnergyFlux_h
#define KM3HAEnergyFlux_h 1

#include <vector>
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "KM3EMDistanceFlux.hh"

class KM3HAEnergyFlux {
public:
  KM3HAEnergyFlux(char *, double, double, double, double);
  ~KM3HAEnergyFlux();

public:
  void FindBins(int idbeam, double energyin, double distancein,
                double anglein);
  int GetNumberOfSamples() { return NumberOfSamples; };
  onePE GetSamplePoint();

private:
  double Rescale(int idbeam, double energyin, double energydis);
  int NPartsDists;
  std::vector<KM3EMDistanceFlux *> *keepEnergies;
  double ParticleEnergyScaleFactors[7][5];
  int ibin;
  double Flux;
  double FluxRMS;
  double EnergyMin;
  double EnergyMax;
  int NumberOfSamples;
  double RatioThis;
};

#endif
