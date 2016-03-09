#ifndef KM3EMShowerModel_h
#define KM3EMShowerModel_h 1

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
#include "KM3EMEnergyFlux.hh"
#include "KM3SD.hh"

class KM3EMShowerModel : public G4VFastSimulationModel {
 public:
  KM3EMShowerModel(G4String, G4Region *);
  KM3EMShowerModel(G4String);
  ~KM3EMShowerModel();

  G4bool IsApplicable(const G4ParticleDefinition &);
  G4bool ModelTrigger(const G4FastTrack &);
  void DoIt(const G4FastTrack &, G4FastStep &);
  KM3Detector *myStDetector;
  KM3SD *aMySD;
  void InitializeFlux(char *, G4double, G4double);

 private:
  KM3EMEnergyFlux *myFlux;
  G4double EnergyThreshold;
  G4double thespeedmaxQE;
};
#endif
