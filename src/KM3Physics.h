#ifndef KM3Physics_h
#define KM3Physics_h 1

#include "G4VUserPhysicsList.h"
#include "globals.h"
#include "KM3Detector.h"
#include "KM3Cherenkov.h"

class KM3Physics : public G4VUserPhysicsList {
 public:
  KM3Physics();
  ~KM3Physics();

 protected:
  // Construct particle and physics
  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();

 protected:
  // these methods Construct physics processes and register them
  void AddParameterisation();
  void ConstructGeneral();
  void ConstructEM();
  void ConstructOP();

#ifdef G4HADRONIC_COMPILE
 protected:
  void ConstructHA();
#endif

 public:
  KM3Detector *aDetector;
  KM3Cherenkov *theCerenkovProcess;

 protected:
  G4double defaultCutEnergyValueForGamma;
  G4double defaultCutEnergyValueForElectron;
  G4double defaultCutEnergyValueForMuon;
  G4double defaultCutEnergyValueForHadron;
};

#endif
