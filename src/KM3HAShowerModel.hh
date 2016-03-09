#ifndef KM3HAShowerModel_h
#define KM3HAShowerModel_h 1

#include "math.h"
#include "G4VFastSimulationModel.hh"
#include "G4Step.hh"
#include "G4TouchableHandle.hh"
#include "G4OpticalPhoton.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "KM3Detector.hh"
#include <vector>
#include "KM3HAEnergyFlux.hh"
#include "KM3SD.hh"

class KM3HAShowerModel : public G4VFastSimulationModel {
 public:
  KM3HAShowerModel(G4String, G4Region *);
  KM3HAShowerModel(G4String);
  ~KM3HAShowerModel();


  G4bool IsApplicable(const G4ParticleDefinition &);
  G4bool ModelTrigger(const G4FastTrack &);
  void DoIt(const G4FastTrack &, G4FastStep &);
  KM3Detector *myStDetector;
  KM3SD *aMySD;
  void InitializeFlux(char *, G4double, G4double);

 private:
  KM3HAEnergyFlux *myFlux;
  G4double EnergyMin;
  G4double EnergyMax;
  G4double thespeedmaxQE;
};
#endif
