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

 private:
  KM3Detector *MyStDetector;

 protected:
};

#endif
