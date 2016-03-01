#ifndef KM3EMDistanceFlux_h
#define KM3EMDistanceFlux_h 1

#include <vector>
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "KM3EMAngularFlux.hh"

class KM3EMDistanceFlux {
public:
  KM3EMDistanceFlux(std::ifstream &);
  ~KM3EMDistanceFlux();

public:
  void FindBins(double distancein, double anglein);
  double GiveEnergy() { return Energy; };
  double GiveFlux() { return Flux; };
  double GiveFluxRMS() { return FluxRMS; };
  onePE GetSamplePoint();

private:
  std::vector<KM3EMAngularFlux *> *keepDistances;
  double Energy;
  int ibin1;
  int ibin2;
  double Flux;
  double FluxRMS;
  double ratio;
  int VertexDistanceBins;
};

#endif
