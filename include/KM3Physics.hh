// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02PhysicsList.hh,v 1.7 2000/12/04 16:24:05 maire Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef KM3Physics_h
#define KM3Physics_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"
#include "KM3Detector.hh"
#include "KM3Cherenkov.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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
