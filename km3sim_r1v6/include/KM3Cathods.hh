#ifndef KM3Cathods_h
#define KM3Cathods_h 1

#include <vector>
#include <stdio.h>
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

struct Cathod {
  G4ThreeVector Position;
  G4ThreeVector Direction;
  G4double Radius;
  G4double Height;
  G4Transform3D trans;
  G4int Depth;
  std::vector<G4int> *Tree;
};

class KM3Cathods {
public:
  KM3Cathods();
  ~KM3Cathods();

public:
  void addCathod(const G4Transform3D &, const G4ThreeVector &,
                 const G4ThreeVector &, const G4double, const G4double,
                 const G4int);
  void addToTree(const G4int);

  G4int GetCathodId(const G4int, const G4int[]);
  void PrintAllCathods(FILE *);
  inline G4Transform3D GetTransformation();
  inline G4ThreeVector GetDirection();
  inline G4ThreeVector GetDirection(G4int it);
  inline G4ThreeVector GetPosition();
  inline G4ThreeVector GetPosition(G4int it);
  inline G4double GetCathodRadius();
  inline G4double GetCathodRadius(G4int it);
  inline G4double GetCathodHeight();
  inline G4double GetCathodHeight(G4int it);
  inline G4int GetNumberOfCathods();

private:
  std::vector<Cathod *> theCathods;
  G4int NumOfCathods;
  G4int iterator;
};

inline G4Transform3D KM3Cathods::GetTransformation() {
  return theCathods[iterator]->trans;
}
inline G4ThreeVector KM3Cathods::GetDirection() {
  return theCathods[iterator]->Direction;
}
inline G4ThreeVector KM3Cathods::GetDirection(G4int it) {
  return theCathods[it]->Direction;
}
inline G4ThreeVector KM3Cathods::GetPosition() {
  return theCathods[iterator]->Position;
}
inline G4ThreeVector KM3Cathods::GetPosition(G4int it) {
  return theCathods[it]->Position;
}
inline G4double KM3Cathods::GetCathodRadius() {
  return theCathods[iterator]->Radius;
}
inline G4double KM3Cathods::GetCathodRadius(G4int it) {
  return theCathods[it]->Radius;
}
inline G4double KM3Cathods::GetCathodHeight() {
  return theCathods[iterator]->Height;
}
inline G4double KM3Cathods::GetCathodHeight(G4int it) {
  return theCathods[it]->Height;
}
inline G4int KM3Cathods::GetNumberOfCathods() { return NumOfCathods; }

#endif
