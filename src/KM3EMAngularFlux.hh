#ifndef KM3EMAngularFlux_h
#define KM3EMAngularFlux_h 1

#include <vector>
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "KM3EMTimePointDis.hh"

class KM3EMAngularFlux {
public:
  KM3EMAngularFlux(std::ifstream &, bool &ok, bool FineBin);
  ~KM3EMAngularFlux();

public:
  void FindBins(double anglein);
  double GiveDistance() { return Distance; };
  double GiveFlux() { return Flux; };
  double GiveFluxRMS() { return FluxRMS; };
  bool IsValid() { return IsThisValid; };
  onePE GetSamplePoint();

private:
  std::vector<KM3EMTimePointDis *> *keepAngles;
  int VertexSolidAngleBins;
  double Distance;
  bool IsThisValid;
  int ibin1;
  int ibin2;
  double Flux;
  double FluxRMS;
  double ratio;
};

#endif
