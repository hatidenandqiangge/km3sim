#ifndef KM3StackingAction_H
#define KM3StackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "KM3Detector.hh"
#include "KM3Cherenkov.hh"
#include "KM3EMDeltaFlux.hh"


class KM3StackingAction : public G4UserStackingAction {
 public:
  KM3StackingAction();
  virtual ~KM3StackingAction();

 public:
  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track *aTrack);
  virtual void NewStage();
  virtual void PrepareNewEvent();
  void SetDetector(KM3Detector *);
#ifdef G4MYHAMUONS_PARAMETERIZATION
  std::ofstream *outMuonHAFile;
#endif

 private:
  KM3Detector *MyStDetector;
#if !defined(G4MYEM_PARAMETERIZATION) && !defined(G4MYHA_PARAMETERIZATION)
#ifdef G4ENABLE_MIE
#ifndef G4DISABLE_PARAMETRIZATION
  // delta rays parametrization section. Works only for muons primaries (not for
  // showers ve)
  std::vector<G4int> *idprikeep;
  std::vector<G4double> *depenekeep;
  std::vector<G4double> *timekeep;
  std::vector<G4ThreeVector> *poskeep;
  // before  G4int indexkeep;
  KM3EMDeltaFlux *myFlux;
  void CreateAllWaitingPhotons(void);
////////////////////////////////////
#endif
#endif
#endif

 protected:
};

#endif
