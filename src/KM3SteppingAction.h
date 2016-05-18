#ifndef KM3SteppingAction_h
#define KM3SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "KM3Detector.h"
#include "KM3EventAction.h"

class KM3SteppingAction : public G4UserSteppingAction {
 public:
  KM3SteppingAction();
  virtual ~KM3SteppingAction(){};
  virtual void UserSteppingAction(const G4Step *aStep);
  KM3Detector *myStDetector;
  KM3EventAction *event_action;

 private:
  G4double MuonRange(G4double);
  G4double P7[8];
  G4double P1LOW[2];
  G4double P1HIGH[2];
};

#endif
