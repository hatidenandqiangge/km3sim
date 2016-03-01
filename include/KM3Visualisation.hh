#ifndef KM3Visualisation_h
#define KM3Visualisation_h 1

#ifdef G4VIS_USE

#include "G4VisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class KM3Visualisation : public G4VisManager {

public:
  KM3Visualisation();

private:
  void RegisterGraphicsSystems();
};

#endif

#endif
