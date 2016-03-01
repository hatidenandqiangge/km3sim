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
  double Radius;
  double Height;
  G4Transform3D trans;
  int Depth;
  std::vector<int> *Tree;
};

class KM3Cathods {
public:
  KM3Cathods();
  ~KM3Cathods();

public:
  void addCathod(const G4Transform3D &, const G4ThreeVector &,
                 const G4ThreeVector &, const double, const double,
                 const int);
  void addToTree(const int);

  int GetCathodId(const int, const int[]);
  void PrintAllCathods(FILE *);
  inline G4Transform3D GetTransformation();
  inline G4ThreeVector GetDirection();
  inline G4ThreeVector GetDirection(int it);
  inline G4ThreeVector GetPosition();
  inline G4ThreeVector GetPosition(int it);
  inline double GetCathodRadius();
  inline double GetCathodRadius(int it);
  inline double GetCathodHeight();
  inline double GetCathodHeight(int it);
  inline int GetNumberOfCathods();

private:
  std::vector<Cathod *> theCathods;
  int NumOfCathods;
  int iterator;
};

inline G4Transform3D KM3Cathods::GetTransformation() {
  return theCathods[iterator]->trans;
}
inline G4ThreeVector KM3Cathods::GetDirection() {
  return theCathods[iterator]->Direction;
}
inline G4ThreeVector KM3Cathods::GetDirection(int it) {
  return theCathods[it]->Direction;
}
inline G4ThreeVector KM3Cathods::GetPosition() {
  return theCathods[iterator]->Position;
}
inline G4ThreeVector KM3Cathods::GetPosition(int it) {
  return theCathods[it]->Position;
}
inline double KM3Cathods::GetCathodRadius() {
  return theCathods[iterator]->Radius;
}
inline double KM3Cathods::GetCathodRadius(int it) {
  return theCathods[it]->Radius;
}
inline double KM3Cathods::GetCathodHeight() {
  return theCathods[iterator]->Height;
}
inline double KM3Cathods::GetCathodHeight(int it) {
  return theCathods[it]->Height;
}
inline int KM3Cathods::GetNumberOfCathods() { return NumOfCathods; }

#endif
