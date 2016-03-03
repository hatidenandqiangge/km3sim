#ifndef KM3MuonParam_h
#define KM3MuonParam_h 1

#include <vector>
#include "globals.hh"
#include "Randomize.hh"

struct PDFSList {
  G4double LogEnergy;
  G4double Distance;
  G4double Prob0;
  G4double StopFirstBin;
  CLHEP::RandGeneral *thePDF;
};

struct forMuon {
  G4double LogEnergy;
  G4double Distance;
  G4double Prob0;
  G4double StopFirstBin;
  CLHEP::RandGeneral *thePDF;
  G4bool iscapable;
};

class KM3MuonParam {
 public:
  KM3MuonParam();
  ~KM3MuonParam();
  G4bool IsCapable(G4int idmuon);
  G4double GetEnergy(G4int idmuon);
  G4double GetDistance(G4int idmuon);
  void AddMuon(G4double energy, G4double distance);
  G4double GetWeight(void);
  void Finalize(void);
  void Initialize(void);

 private:
  std::vector<PDFSList *> thePDFS;
  std::vector<forMuon *> theDistributions;
  G4double MinLogEnergy;
  G4double MaxLogEnergy;
  G4double MyTolerance;
  G4double EnergyCutoff;
  G4bool isEventOK;
};
#endif
