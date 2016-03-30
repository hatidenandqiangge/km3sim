#ifndef KM3TrackingAction_h
#define KM3TrackingAction_h 1

#include "G4UserTrackingAction.h"
#include "HOURSevtWRITE.h"
#include "G4Types.h"

class KM3TrackingAction : public G4UserTrackingAction {
 public:
  KM3TrackingAction(){};
  virtual ~KM3TrackingAction(){};

  virtual void PreUserTrackingAction(const G4Track *);
  virtual void PostUserTrackingAction(const G4Track *);

 public:
  int numofInitialParticles;
  HOURSevtWRITE *TheEVTtoWrite;
  G4bool useANTARESformat;
};

#endif
