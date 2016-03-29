#ifndef KM3EventAction_h
#define KM3EventAction_h 1
#include "stdio.h"
#include <vector>
#include "G4Types.hh"
#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"

#include "HOURSevtWRITE.hh"

class G4EventManager;
class G4Event;

// class description:
//
//  This is the base class of one of the user's optional action classes.
// The two methods BeginOfEventAction() and EndOfEventAction() are invoked
// at the beginning and the end of one event processing. These methods are
// invoked by G4EventManager.
//  Be aware that BeginOfEventAction() is invoked when a G4Event object is
// sent to G4EventManager. Thus the primary vertexes/particles have already
// been made by the primary generator. In case the user wants to do something
// before generating primaries (i.e., store random number status), do it in
// the G4VUserPrimaryGeneratorAction concrete class.
//

class KM3EventAction : public G4UserEventAction {
 public:
  KM3EventAction() { ; }
  ~KM3EventAction() { ; }
  inline void SetEventManager(G4EventManager *value) { fpEventManager = value; }

 public:  // with description
  void BeginOfEventAction(const G4Event *anEvent);
  void EndOfEventAction(const G4Event *anEvent);
  // Two virtual method the user can override.
 protected:
  G4EventManager *fpEventManager;

 public:
  std::vector<G4ThreeVector> centerPre;
  std::vector<G4ThreeVector> centerPost;
  std::vector<G4ThreeVector> enterPre;
  std::vector<G4ThreeVector> enterPost;
  std::vector<G4ThreeVector> leavePre;
  std::vector<G4ThreeVector> leavePost;
  std::vector<G4double> centerMomentum;
  std::vector<G4double> enterMomentum;
  std::vector<G4double> leaveMomentum;
  std::vector<G4ThreeVector> centerPosition;
  std::vector<G4ThreeVector> enterPosition;
  std::vector<G4ThreeVector> leavePosition;
  std::vector<G4double> centerTime;
  std::vector<G4double> enterTime;
  std::vector<G4double> leaveTime;
  std::vector<G4ThreeVector> stopPosition;
  std::vector<G4double> stopTime;
#ifdef G4MYMUON_KEEPENERGY
  std::vector<G4double> EnergyAtPosition;
#endif
  G4bool useANTARESformat;
  HOURSevtWRITE *TheEVTtoWrite;

 public:
  FILE *outfile;
  inline void AddPrimaryNumber(G4int);
  inline G4int GetSlot(G4int);
  inline void Initialize(void) { numofMuons = 0; }

 private:
  G4int numofMuons;
  G4int MuonIds[210000];  // 100 before
};

inline void KM3EventAction::AddPrimaryNumber(G4int aNumber) {
  MuonIds[numofMuons] = aNumber;
  numofMuons++;
}

inline G4int KM3EventAction::GetSlot(G4int aNumber) {
  G4int slot = -1;
  for (G4int i = 0; i < numofMuons; i++) {
    if (MuonIds[i] == aNumber) {
      slot = i;
      break;
    }
  }
  return slot;
}

#endif
