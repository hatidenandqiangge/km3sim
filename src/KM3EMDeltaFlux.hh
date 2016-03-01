#ifndef KM3EMDeltaFlux_h
#define KM3EMDeltaFlux_h 1

#include <vector>
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "KM3EMAngularFlux.hh"

class KM3EMDeltaFlux {
public:
  KM3EMDeltaFlux(char *, double, double);
  ~KM3EMDeltaFlux();

public:
  void FindBins(double MeanNumPhotons, double distancein, double anglein);
  int GetNumberOfSamples() { return NumberOfSamples; };
  onePE GetSamplePoint();

private:
  std::vector<KM3EMAngularFlux *> *keepDistances;
  int ibin1;
  int ibin2;
  double Flux;
  double ratio;
  int NumberOfSamples;
  double RatioThis;
  int VertexDistanceBins;
};

#endif
