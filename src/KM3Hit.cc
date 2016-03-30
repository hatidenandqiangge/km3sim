#include "KM3Hit.h"
#include "G4UnitsTable.h"
#include "G4VisManager.h"
#include "G4Circle.h"
#include "G4Colour.h"
#include "G4VisAttributes.h"
#include "G4Transform3D.h"
#include "G4LogicalVolume.h"

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
