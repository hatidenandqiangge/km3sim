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
  KM3EMDirectFlux(char *, double);
  ~KM3EMDirectFlux();

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
