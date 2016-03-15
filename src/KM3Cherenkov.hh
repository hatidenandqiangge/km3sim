#ifndef KM3Cherenkov_H
#define KM3Cherenkov_H 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4PhysicsTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "KM3Detector.hh"
// Class Description:
// Continuous Process -- Generation of Cerenkov Photons.
// Class inherits publicly from G4VDiscreteProcess.
// Class Description - End:

class KM3Cherenkov : public G4VProcess {
 public:
  KM3Cherenkov(const G4String &processName = "KM3Cherenkov",
               G4ProcessType type = fElectromagnetic);

  KM3Cherenkov(const KM3Cherenkov &right);
  ~KM3Cherenkov();

 private:

  KM3Cherenkov &operator=(const KM3Cherenkov &right);

 public:  // With description
  void SetDetector(KM3Detector *);

  G4bool IsApplicable(const G4ParticleDefinition &aParticleType);
  // Returns true -> 'is applicable', for all charged particles.
  // except short-lived particles.

  G4double GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition *);
  // Returns the discrete step limit and sets the 'StronglyForced'
  // condition for the DoIt to be invoked at every step.

  G4double PostStepGetPhysicalInteractionLength(const G4Track &aTrack, G4double,
                                                G4ForceCondition *);
  // Returns the discrete step limit and sets the 'StronglyForced'
  // condition for the DoIt to be invoked at every step.

  G4VParticleChange *PostStepDoIt(const G4Track &aTrack, const G4Step &aStep);
  // This is the method implementing the Cerenkov process.

  //  no operation in  AtRestDoIt and  AlongStepDoIt
  virtual G4double AlongStepGetPhysicalInteractionLength(const G4Track &,
                                                         G4double, G4double,
                                                         G4double &,
                                                         G4GPILSelection *) {
    return -1.0;
  };

  virtual G4double AtRestGetPhysicalInteractionLength(const G4Track &,
                                                      G4ForceCondition *) {
    return -1.0;
  };

  //  no operation in  AtRestDoIt and  AlongStepDoIt
  virtual G4VParticleChange *AtRestDoIt(const G4Track &, const G4Step &) {
    return 0;
  };

  virtual G4VParticleChange *AlongStepDoIt(const G4Track &, const G4Step &) {
    return 0;
  };

  void SetTrackSecondariesFirst(const G4bool state);
  // If set, the primary particle tracking is interrupted and any
  // produced Cerenkov photons are tracked next. When all have
  // been tracked, the tracking of the primary resumes.

  void SetMaxBetaChangePerStep(const G4double d);
  // Set the maximum allowed change in beta = v/c in % (perCent)
  // per step.

  void SetMaxNumPhotonsPerStep(const G4int NumPhotons);
  // Set the maximum number of Cerenkov photons allowed to be
  // generated during a tracking step. This is an average ONLY;
  // the actual number will vary around this average. If invoked,
  // the maximum photon stack will roughly be of the size set.
  // If not called, the step is not limited by the number of
  // photons generated.

  G4PhysicsTable *GetPhysicsTable() const;
  // Returns the address of the physics table.

  void DumpPhysicsTable() const;
// Prints the physics table.

 private:

#ifdef G4JUST_COUNT_PHOTONS
  long double Count_Photons;
  long double Posit_Photons_Mean;
  long double
      Posit_Photons_Histo[3002];  // 20 meters every 1cm from -10m to 20m
                                  // , and 2 for underflow and overflow
#endif

  void BuildThePhysicsTable();

  G4double GetAverageNumberOfPhotons(const G4double charge, const G4double beta,
                                     const G4Material *aMaterial,
                                     G4MaterialPropertyVector *Rindex) const;

 protected:
  G4PhysicsTable *thePhysicsTable;
  //  A Physics Table can be either a cross-sections table or
  //  an energy table (or can be used for other specific
  //  purposes).

 private:
  G4bool fTrackSecondariesFirst;
  G4double fMaxBetaChange;
  G4int fMaxPhotons;
  KM3Detector *MyStDetector;
  G4double MaxAbsDist;
  G4double M_PI2;
  G4double MinMeanNumberOfPhotonsForParam;
};

inline G4bool KM3Cherenkov::IsApplicable(
    const G4ParticleDefinition &aParticleType) {
  if (aParticleType.GetParticleName() == "chargedgeantino") return false;
  if (aParticleType.IsShortLived()) return false;
  return (aParticleType.GetPDGCharge() != 0);
}

inline void KM3Cherenkov::SetTrackSecondariesFirst(const G4bool state) {
  fTrackSecondariesFirst = state;
}

inline void KM3Cherenkov::SetMaxBetaChangePerStep(const G4double value) {
  fMaxBetaChange = value * CLHEP::perCent;
}

inline void KM3Cherenkov::SetMaxNumPhotonsPerStep(const G4int NumPhotons) {
  fMaxPhotons = NumPhotons;
}

inline void KM3Cherenkov::DumpPhysicsTable() const {
  G4int PhysicsTableSize = thePhysicsTable->entries();
  G4PhysicsOrderedFreeVector *v;

  for (G4int i = 0; i < PhysicsTableSize; i++) {
    v = (G4PhysicsOrderedFreeVector *)(*thePhysicsTable)[i];
    v->DumpValues();
  }
}

inline G4PhysicsTable *KM3Cherenkov::GetPhysicsTable() const {
  return thePhysicsTable;
}

#endif /* G4Cerenkov_h */
