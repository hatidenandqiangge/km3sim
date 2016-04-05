#ifndef KM3PrimaryGeneratorAction_h
#define KM3PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.h"
#include "KM3TrackingAction.h"
#include "G4ThreeVector.h"
#include "KM3EventAction.h"

#include <stdio.h>
#include "globals.h"

#include "HAVertexMuons.h"

class G4Event;
class G4VPrimaryGenerator;

#include "HOURSevtREAD.h"

class KM3PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  KM3PrimaryGeneratorAction();
  ~KM3PrimaryGeneratorAction();

 public:
  int nevents;
  FILE *outfile;
  FILE *infile;
  char *filePythiaParticles;
  char *fileParticles;
  G4int numberofParticles;
  void GeneratePrimaries(G4Event *anEvent);
  void Initialize(void);
  G4bool useHEPEvt;
  G4bool useANTARESformat;
  KM3TrackingAction *myTracking;
  KM3EventAction *event_action;
  G4double ParamEnergy;
  G4int idbeam;  // type of injected or produced particles (PDG Code)
  G4double random_R;
  G4ThreeVector position;
  G4ThreeVector direction;

 private:
  G4VPrimaryGenerator *HEPEvt;
  HOURSevtREAD *antaresHEPEvt;
  G4double EventWeight;
  G4ThreeVector detectorCenter;
  G4double detectorMaxRho;
  G4double detectorMaxz;
  G4double bottomPosition;

  HAVertexMuons *aHAVertexMuons;

 public:
  void PutFromDetector(G4ThreeVector dC, G4double dMR, G4double dMz,
                       G4double bP) {
    detectorCenter = dC;
    detectorMaxRho = dMR;
    detectorMaxz = dMz;
    bottomPosition = bP;
  };
};

#endif /*KM3PrimaryGeneratorAction_h*/
