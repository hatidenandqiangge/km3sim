//
#include "KM3Hit.hh"
#include "G4UnitsTable.hh"
#include "G4VisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"

G4Allocator<KM3Hit> KM3HitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KM3Hit::KM3Hit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KM3Hit::~KM3Hit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KM3Hit::KM3Hit(const KM3Hit& right)
{
  CathodId  = right.CathodId;
  time      = right.time;
  originalInfo=right.originalInfo;
  IMany=right.IMany;
  //short  angleIncident=right.angleIncident;
  //short  angleDirection=right.angleDirection;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const KM3Hit& KM3Hit::operator=(const KM3Hit& right)
{
  CathodId  = right.CathodId;
  time      = right.time;
  originalInfo=right.originalInfo;
  IMany=right.IMany;
  //short  angleIncident=right.angleIncident;
  //short  angleDirection=right.angleDirection;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int KM3Hit::operator==(const KM3Hit& right) const
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

