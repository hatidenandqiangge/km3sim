#include "KM3Hit.h"
#include "G4UnitsTable.hh"
#include "G4VisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"

G4Allocator<KM3Hit> KM3HitAllocator;

KM3Hit::KM3Hit() {}

KM3Hit::~KM3Hit() {}

KM3Hit::KM3Hit(const KM3Hit &right) {
  CathodId = right.CathodId;
  time = right.time;
  originalInfo = right.originalInfo;
  IMany = right.IMany;
  // short  angleIncident=right.angleIncident;
  // short  angleDirection=right.angleDirection;
}

const KM3Hit &KM3Hit::operator=(const KM3Hit &right) {
  CathodId = right.CathodId;
  time = right.time;
  originalInfo = right.originalInfo;
  IMany = right.IMany;
  // short  angleIncident=right.angleIncident;
  // short  angleDirection=right.angleDirection;
  return *this;
}

int KM3Hit::operator==(const KM3Hit &right) const { return 0; }
