#ifndef KM3MuonParam_h
#define KM3MuonParam_h 1

#include <vector>
#include "globals.hh"
#include "Randomize.hh"

struct PDFSList {
  double LogEnergy;
  double Distance;
  double Prob0;
  double StopFirstBin;
  CLHEP::RandGeneral *thePDF;
};

struct forMuon {
  double LogEnergy;
  double Distance;
  double Prob0;
  double StopFirstBin;
  CLHEP::RandGeneral *thePDF;
  bool iscapable;
};

class KM3MuonParam {
public:
  KM3MuonParam();
  ~KM3MuonParam();
  bool IsCapable(int idmuon);
  double GetEnergy(int idmuon);
  double GetDistance(int idmuon);
  void AddMuon(double energy, double distance);
  double GetWeight(void);
  void Finalize(void);
  void Initialize(void);

private:
  std::vector<PDFSList *> thePDFS;
  std::vector<forMuon *> theDistributions;
  double MinLogEnergy;
  double MaxLogEnergy;
  double MyTolerance;
  double EnergyCutoff;
  bool isEventOK;
};
#endif
