#ifndef HAVertexMuons_h
#define HAVertexMuons_h 1

#include <vector>
#include "globals.hh"
#include "G4ThreeVector.hh"
#include <iostream>
#include <fstream>
#include <iomanip>

class HAVertexMuons {
public:
  HAVertexMuons(char *, char *);
  ~HAVertexMuons();
  int GetNumberOfMuons(double HadronicEnergy);
  void ReadMuon();
  inline G4ThreeVector GetPosition() { return thePosition; };
  inline G4ThreeVector GetMomentum() { return theMomentum; };
  inline double GetTime() { return theTime; };

private:
  std::ifstream *MuonsStream;
  int numevents;
  std::vector<float> *theEnergies;
  std::vector<int> *theSizes;
  std::vector<int> *thePositions;
  double theTime;
  G4ThreeVector thePosition;
  G4ThreeVector theMomentum;
};

#endif
