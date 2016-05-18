#ifndef HAVertexMuons_h
#define HAVertexMuons_h 1

#include "globals.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "G4ThreeVector.hh"

class HAVertexMuons {
 public:
  HAVertexMuons(char *, char *);
  ~HAVertexMuons();
  G4int GetNumberOfMuons(G4double HadronicEnergy);
  void ReadMuon();
  inline G4ThreeVector GetPosition() { return thePosition; };
  inline G4ThreeVector GetMomentum() { return theMomentum; };
  inline G4double GetTime() { return theTime; };

 private:
  std::ifstream *MuonsStream;
  G4int numevents;
  std::vector<G4float> *theEnergies;
  std::vector<G4int> *theSizes;
  std::vector<G4int> *thePositions;
  G4double theTime;
  G4ThreeVector thePosition;
  G4ThreeVector theMomentum;
};

#endif
