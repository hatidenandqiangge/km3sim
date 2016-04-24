#include "G4ios.h"
#include "G4PhysicalConstants.h"
#include "G4SystemOfUnits.h"
#include "G4Poisson.h"
#include "G4EmProcessSubType.h"
#include "G4LossTableManager.h"
#include "G4MaterialCutsCouple.h"
#include "G4ParticleDefinition.h"
#include "KM3Cherenkov.h"

KM3Cherenkov::KM3Cherenkov(const G4String &processName, G4ProcessType type)
    : G4VProcess(processName, type) {
  SetProcessSubType(fCerenkov);

  fTrackSecondariesFirst = false;
  fMaxBetaChange = 0.;
  fMaxPhotons = 0;

  thePhysicsTable = NULL;

  if (verboseLevel > 0) {
    G4cout << GetProcessName() << " is created " << G4endl;
  }
  BuildThePhysicsTable();

  M_PI2 = 2 * M_PI;
  MinMeanNumberOfPhotonsForParam = 20.0;

#ifdef G4JUST_COUNT_PHOTONS
  Count_Photons = 0.0;
  Posit_Photons_Mean = 0.0;
  for (G4int i = 0; i < 3002; i++) Posit_Photons_Histo[i] = 0.0;
#endif
}

KM3Cherenkov::~KM3Cherenkov() {
  if (thePhysicsTable != NULL) {
    thePhysicsTable->clearAndDestroy();
    delete thePhysicsTable;
  }

#ifdef G4JUST_COUNT_PHOTONS
  G4int ibin;
  long double cumul = 0;
  for (ibin = 0; ibin < 3002; ibin++) {
    cumul += Posit_Photons_Histo[ibin];
    if (cumul > 0.5 * Count_Photons) break;
  }
  G4double Posit_Photons_Median = -10 * m + 0.5 * cm + G4double(ibin - 1) * cm;
  Posit_Photons_Mean /= Count_Photons;
  printf("Count_Photons %.20Le %.5Le %.5e\n", Count_Photons,
         Posit_Photons_Mean / m, Posit_Photons_Median / m);
#endif
}

void KM3Cherenkov::SetDetector(KM3Detector *adet) {
  MyStDetector = adet;
  MaxAbsDist = MyStDetector->MaxAbsDist;
}

// This is the method implementing the Cerenkov process.
G4VParticleChange *KM3Cherenkov::PostStepDoIt(const G4Track &aTrack,
                                              const G4Step &aStep) {
  // G4cout << " IN CHERENKOV";
  // G4cout << " Particle id " << aTrack.GetTrackID() << G4endl;
  // G4cout << " Parent id " << aTrack.GetParentID();
  // G4cout << " Particle name " << aTrack.GetDefinition()->GetParticleName();
  // G4cout << " Kinetic Energy " << aTrack.GetKineticEnergy() << " " <<
  // aTrack.GetDynamicParticle()->GetKineticEnergy();
  // G4cout << " Total Energy " << aTrack.GetTotalEnergy();
  // G4cout << " Creator Process " <<
  // aTrack.GetCreatorProcess()->GetProcessName();
  // G4cout << " Track Length " << aTrack.GetTrackLength();
  // G4cout << " STEP Length from step " << aStep.GetStepLength();
  // G4cout << " delta position phi " << aStep.GetDeltaPosition().phi()*deg;
  // G4cout << " delta position theta " << aStep.GetDeltaPosition().theta()*deg;
  // G4cout << " delta position magnitude " <<
  // aStep.GetDeltaPosition().mag()<<G4endl;

  /// at first initialize the pointers to Q_E, glass and gell transparencies////
  static G4MaterialPropertyVector *QECathod = NULL;
  if (QECathod == NULL) {
    const G4MaterialTable *theMaterialTable = G4Material::GetMaterialTable();
    for (size_t J = 0; J < theMaterialTable->size(); J++) {
      if ((*theMaterialTable)[J]->GetName() == G4String("Cathod")) {
        G4MaterialPropertiesTable *aMaterialPropertiesTable =
            (*theMaterialTable)[J]->GetMaterialPropertiesTable();
        QECathod = aMaterialPropertiesTable->GetProperty("Q_EFF");
      }
    }
  }

  // Should we ensure that the material is dispersive?
  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle *aParticle = aTrack.GetDynamicParticle();
  const G4Material *aMaterial = aTrack.GetMaterial();

  G4MaterialPropertiesTable *aMaterialPropertiesTable =
      aMaterial->GetMaterialPropertiesTable();

  if (!aMaterialPropertiesTable) return pParticleChange;

  G4MaterialPropertyVector *Rindex =
      aMaterialPropertiesTable->GetProperty("RINDEX");

  if (!Rindex) return pParticleChange;

  // check that the particle is inside the active volume of the detector
  static G4double detectorMaxRho2 =
      MyStDetector->detectorMaxRho * MyStDetector->detectorMaxRho;
  G4StepPoint *pPreStepPoint = aStep.GetPreStepPoint();
  G4ThreeVector x0 = pPreStepPoint->GetPosition();
  //  G4cout <<"prepoint "<< x0[0] <<" "<< x0[1] <<" "<< x0[2] <<G4endl;
  G4ThreeVector distanceV = x0 - MyStDetector->detectorCenter;
  G4double distanceRho2 =
      distanceV[0] * distanceV[0] + distanceV[1] * distanceV[1];
  if ((distanceRho2 > detectorMaxRho2) ||
      (x0[2] < MyStDetector->bottomPosition) ||
      (x0[2] > MyStDetector->detectorMaxz)) {
    // return unchanged particle and no secondaries
    aParticleChange.SetNumberOfSecondaries(0);
    return pParticleChange;
  }

  // step length
  G4double step_length = aStep.GetStepLength();
  //  G4cout << "step_length "<<aParticle->GetDefinition()->GetParticleName()
  //  <<" "<< step_length<<G4endl;
  // particle charge
  const G4double charge = aParticle->GetDefinition()->GetPDGCharge();
  // particle beta
  G4StepPoint *pPostStepPoint = aStep.GetPostStepPoint();
  const G4double beta =
      (pPreStepPoint->GetBeta() + pPostStepPoint->GetBeta()) / 2.;
  G4double MeanNumPhotons =
      GetAverageNumberOfPhotons(charge, beta, aMaterial, Rindex);
  MeanNumPhotons *= step_length * MyStDetector->Quantum_Efficiency;
  G4double BetaInverse = 1.0 / beta;
  G4ThreeVector p0 = aStep.GetDeltaPosition().unit();

  G4bool EmittedAsScattered = false;  // newmie

  if (MeanNumPhotons <= 0.0) {
    // return unchanged particle and no secondaries
    aParticleChange.SetNumberOfSecondaries(0);
    return pParticleChange;
  }

  G4int NumPhotons = (G4int)CLHEP::RandPoisson::shoot(MeanNumPhotons);
  if (NumPhotons <= 0) {
    // return unchanged particle and no secondaries
    aParticleChange.SetNumberOfSecondaries(0);
    return pParticleChange;
  }

#ifdef G4JUST_COUNT_PHOTONS
  Count_Photons += (long double)NumPhotons;
  G4double pos_z_projection =
      0.5 * ((aStep.GetPreStepPoint()->GetPosition())[2] +
             (aStep.GetPostStepPoint()->GetPosition())[2]);
  Posit_Photons_Mean +=
      ((long double)NumPhotons) * ((long double)pos_z_projection);
  G4int ibin = 1001 + int(floor(pos_z_projection / cm));
  if (ibin < 0) ibin = 0;
  if (ibin > 3001) ibin = 3001;
  Posit_Photons_Histo[ibin] += (long double)NumPhotons;
  aParticleChange.SetNumberOfSecondaries(0);
  return pParticleChange;
#endif

  G4double nMax = Rindex->GetMaxValue();
  G4double maxCos = BetaInverse / nMax;

  aParticleChange.SetNumberOfSecondaries(NumPhotons);

  if (fTrackSecondariesFirst) {
    if (aTrack.GetTrackStatus() == fAlive)
      aParticleChange.ProposeTrackStatus(fSuspend);
  }

  G4double Pmin = Rindex->GetMinLowEdgeEnergy();
  G4double Pmax = Rindex->GetMaxLowEdgeEnergy();
  G4double dp = Pmax - Pmin;
  G4double maxSin2 = (1.0 - maxCos) * (1.0 + maxCos);
  const G4double beta1 = pPreStepPoint->GetBeta();
  const G4double beta2 = pPostStepPoint->GetBeta();
  G4double MeanNumberOfPhotons1 =
      GetAverageNumberOfPhotons(charge, beta1, aMaterial, Rindex);
  G4double MeanNumberOfPhotons2 =
      GetAverageNumberOfPhotons(charge, beta2, aMaterial, Rindex);
  G4double t0 = pPreStepPoint->GetGlobalTime();
  //  NumPhotons=0; //lookout
  for (G4int i = 0; i < NumPhotons; i++) {
    G4double rand;
    G4double sampledEnergy, sampledRI;
    G4double cosTheta, sin2Theta;

    // sample a phi
    rand = G4UniformRand();
    G4double phi = M_PI2 * rand;
    G4double sinPhi = sin(phi);
    G4double cosPhi = cos(phi);

    // Determine photon energy
    // sample an energy

    do {
      rand = G4UniformRand();
      sampledEnergy = Pmin + rand * dp;
      sampledRI = Rindex->Value(sampledEnergy);
      cosTheta = BetaInverse / sampledRI;

      sin2Theta = (1.0 - cosTheta) * (1.0 + cosTheta);
      rand = G4UniformRand();

    } while (rand * maxSin2 > sin2Theta);

    G4double qeProb = QECathod->Value(sampledEnergy);

    if (G4UniformRand() < qeProb) {
      // calculate x,y, and z components of photon momentum
      // (in coord system with primary particle direction
      //  aligned with the z axis)
      G4double sinTheta = sqrt(sin2Theta);
      G4double px = sinTheta * cosPhi;
      G4double py = sinTheta * sinPhi;
      G4double pz = cosTheta;

      // Create photon momentum direction vector
      // The momentum direction is still with respect
      // to the coordinate system where the primary
      // particle direction is aligned with the z axis
      G4ParticleMomentum photonMomentum(px, py, pz);

      // Rotate momentum direction back to global reference
      // system
      photonMomentum.rotateUz(p0);

      // Determine polarization of new photon
      G4double sx = cosTheta * cosPhi;
      G4double sy = cosTheta * sinPhi;
      G4double sz = -sinTheta;
      G4ThreeVector photonPolarization(sx, sy, sz);

      // Rotate back to original coord system
      photonPolarization.rotateUz(p0);

      // Generate new G4Track object and a new photon
      // first find generation position and time
      G4double delta, NumberOfPhotons, N;
      do {
        rand = G4UniformRand();
        delta = rand * step_length;
        NumberOfPhotons =
            MeanNumberOfPhotons1 -
            delta * (MeanNumberOfPhotons1 - MeanNumberOfPhotons2) / step_length;
        N = G4UniformRand() *
            std::max(MeanNumberOfPhotons1, MeanNumberOfPhotons2);
      } while (N > NumberOfPhotons);

      G4double deltaTime =
          delta /
          ((pPreStepPoint->GetVelocity() + pPostStepPoint->GetVelocity()) / 2.);

      G4ThreeVector aSecondaryPosition = x0 + rand * aStep.GetDeltaPosition();

      G4double aSecondaryTime = t0 + deltaTime;

      G4DynamicParticle *aCerenkovPhoton = new G4DynamicParticle(
          G4OpticalPhoton::OpticalPhoton(), photonMomentum);
      aCerenkovPhoton->SetPolarization(photonPolarization.x(),
                                       photonPolarization.y(),
                                       photonPolarization.z());

      aCerenkovPhoton->SetKineticEnergy(sampledEnergy);

      // Generate the track
      G4Track *aSecondaryTrack =
          new G4Track(aCerenkovPhoton, aSecondaryTime, aSecondaryPosition);

      aSecondaryTrack->SetTouchableHandle(
          aStep.GetPreStepPoint()->GetTouchableHandle());
      aSecondaryTrack->SetParentID(aTrack.GetTrackID());

      aParticleChange.AddSecondary(aSecondaryTrack);
    }  // if (G4UniformRand()<qeProb)
  }    // for each photon

  if (verboseLevel > 0) {
    G4cout << "\n Exiting from KM3Cherenkov::DoIt -- NumberOfSecondaries = "
           << aParticleChange.GetNumberOfSecondaries() << G4endl;
  }

  return pParticleChange;
}

void KM3Cherenkov::BuildThePhysicsTable() {
  if (thePhysicsTable) return;

  const G4MaterialTable *theMaterialTable = G4Material::GetMaterialTable();
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();

  // create new physics table

  thePhysicsTable = new G4PhysicsTable(numOfMaterials);

  // loop for materials

  for (G4int i = 0; i < numOfMaterials; i++) {
    G4PhysicsOrderedFreeVector *aPhysicsOrderedFreeVector =
        new G4PhysicsOrderedFreeVector();

    // Retrieve vector of refraction indices for the material
    // from the material's optical properties table

    G4Material *aMaterial = (*theMaterialTable)[i];

    G4MaterialPropertiesTable *aMaterialPropertiesTable =
        aMaterial->GetMaterialPropertiesTable();

    if (aMaterialPropertiesTable) {
      G4MaterialPropertyVector *theRefractionIndexVector =
          aMaterialPropertiesTable->GetProperty("RINDEX");

      if (theRefractionIndexVector) {
        // Retrieve the first refraction index in vector
        // of (photon momentum, refraction index) pairs

        G4double currentRI = (*theRefractionIndexVector)[0];

        if (currentRI > 1.0) {
          // Create first (photon momentum, Cerenkov Integral)
          // pair

          G4double currentPM = theRefractionIndexVector->Energy(0);
          G4double currentCAI = 0.0;

          aPhysicsOrderedFreeVector->InsertValues(currentPM, currentCAI);

          // Set previous values to current ones prior to loop

          G4double prevPM = currentPM;
          G4double prevCAI = currentCAI;
          G4double prevRI = currentRI;

          // loop over all (photon momentum, refraction index)
          // pairs stored for this material

          for (size_t i = 1; i < theRefractionIndexVector->GetVectorLength();
               i++) {
            currentRI = (*theRefractionIndexVector)[i];
            currentPM = theRefractionIndexVector->Energy(i);

            currentCAI =
                0.5 * (1.0 / (prevRI * prevRI) + 1.0 / (currentRI * currentRI));

            currentCAI = prevCAI + (currentPM - prevPM) * currentCAI;

            aPhysicsOrderedFreeVector->InsertValues(currentPM, currentCAI);

            prevPM = currentPM;
            prevCAI = currentCAI;
            prevRI = currentRI;
          }
        }
      }
    }

    // The Cerenkov integral for a given material
    // will be inserted in thePhysicsTable
    // according to the position of the material in
    // the material table.

    thePhysicsTable->insertAt(i, aPhysicsOrderedFreeVector);
  }
}

// GetMeanFreePath
// ---------------
//

G4double KM3Cherenkov::GetMeanFreePath(const G4Track &, G4double,
                                       G4ForceCondition *) {
  return 1.;
}

G4double KM3Cherenkov::PostStepGetPhysicalInteractionLength(
    const G4Track &aTrack, G4double, G4ForceCondition *condition) {
  *condition = NotForced;
  G4double StepLimit = DBL_MAX;

  const G4DynamicParticle *aParticle = aTrack.GetDynamicParticle();
  const G4Material *aMaterial = aTrack.GetMaterial();
  const G4MaterialCutsCouple *couple = aTrack.GetMaterialCutsCouple();

  G4double kineticEnergy = aParticle->GetKineticEnergy();
  const G4ParticleDefinition *particleType = aParticle->GetDefinition();
  G4double mass = particleType->GetPDGMass();

  // particle beta
  G4double beta = aParticle->GetTotalMomentum() / aParticle->GetTotalEnergy();
  // particle gamma
  G4double gamma = aParticle->GetTotalEnergy() / mass;

  G4MaterialPropertiesTable *aMaterialPropertiesTable =
      aMaterial->GetMaterialPropertiesTable();

  G4MaterialPropertyVector *Rindex = NULL;

  if (aMaterialPropertiesTable)
    Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");

  G4double nMax;
  if (Rindex) {
    nMax = Rindex->GetMaxValue();
  } else {
    return StepLimit;
  }

  G4double BetaMin = 1. / nMax;
  if (BetaMin >= 1.) return StepLimit;

  G4double GammaMin = 1. / std::sqrt(1. - BetaMin * BetaMin);

  if (gamma < GammaMin) return StepLimit;

  G4double kinEmin = mass * (GammaMin - 1.);

  G4double RangeMin =
      G4LossTableManager::Instance()->GetRange(particleType, kinEmin, couple);
  G4double Range = G4LossTableManager::Instance()->GetRange(
      particleType, kineticEnergy, couple);

  G4double Step = Range - RangeMin;
  if (Step < 1. * um) return StepLimit;

  if (Step > 0. && Step < StepLimit) StepLimit = Step;

  // If user has defined an average maximum number of photons to
  // be generated in a Step, then calculate the Step length for
  // that number of photons.

  if (fMaxPhotons > 0) {
    // particle charge
    const G4double charge = aParticle->GetDefinition()->GetPDGCharge();

    G4double MeanNumberOfPhotons =
        GetAverageNumberOfPhotons(charge, beta, aMaterial, Rindex);

    Step = 0.;
    if (MeanNumberOfPhotons > 0.0) Step = fMaxPhotons / MeanNumberOfPhotons;

    if (Step > 0. && Step < StepLimit) StepLimit = Step;
  }

  // If user has defined an maximum allowed change in beta per step
  if (fMaxBetaChange > 0.) {
    G4double dedx = G4LossTableManager::Instance()->GetDEDX(
        particleType, kineticEnergy, couple);

    G4double deltaGamma = gamma -
                          1. / std::sqrt(1. -
                                         beta * beta * (1. - fMaxBetaChange) *
                                             (1. - fMaxBetaChange));

    Step = mass * deltaGamma / dedx;

    if (Step > 0. && Step < StepLimit) StepLimit = Step;
  }

  *condition = StronglyForced;
  return StepLimit;
}

// GetAverageNumberOfPhotons
// -------------------------
// This routine computes the number of Cerenkov photons produced per
// GEANT-unit (millimeter) in the current medium.
//             ^^^^^^^^^^

G4double KM3Cherenkov::GetAverageNumberOfPhotons(
    const G4double charge, const G4double beta, const G4Material *aMaterial,
    G4MaterialPropertyVector *Rindex) const {
  const G4double Rfact = 369.81 / (eV * cm);

  if (beta <= 0.0) return 0.0;

  G4double BetaInverse = 1. / beta;

  // Vectors used in computation of Cerenkov Angle Integral:
  //   - Refraction Indices for the current material
  //  - new G4PhysicsOrderedFreeVector allocated to hold CAI's
  G4int materialIndex = aMaterial->GetIndex();

  // Retrieve the Cerenkov Angle Integrals for this material
  G4PhysicsOrderedFreeVector *CerenkovAngleIntegrals =
      (G4PhysicsOrderedFreeVector *)((*thePhysicsTable)(materialIndex));

  if (!(CerenkovAngleIntegrals->IsFilledVectorExist())) return 0.0;

  // Min and Max photon momenta
  G4double Pmin = Rindex->GetMinLowEdgeEnergy();
  G4double Pmax = Rindex->GetMaxLowEdgeEnergy();

  // Min and Max Refraction Indices
  G4double nMin = Rindex->GetMinValue();
  G4double nMax = Rindex->GetMaxValue();

  // Max Cerenkov Angle Integral
  G4double CAImax = CerenkovAngleIntegrals->GetMaxValue();

  G4double dp, ge;

  // If n(Pmax) < 1/Beta -- no photons generated

  if (nMax < BetaInverse) {
    //    G4cout<<aParticle->GetTotalEnergy()<<"
    //"<<aParticle->GetKineticEnergy()<<G4endl;
    //    dp = 0;
    //  ge = 0;
    return 0.0;
  }

  // otherwise if n(Pmin) >= 1/Beta -- photons generated

  else if (nMin > BetaInverse) {
    dp = Pmax - Pmin;
    ge = CAImax;
  }

  // If n(Pmin) < 1/Beta, and n(Pmax) >= 1/Beta, then
  // we need to find a P such that the value of n(P) == 1/Beta.
  // Interpolation is performed by the GetPhotonMomentum() and
  // GetProperty() methods of the G4MaterialPropertiesTable and
  // the GetValue() method of G4PhysicsVector.

  else {
    Pmin = Rindex->GetEnergy(BetaInverse);
    dp = Pmax - Pmin;

    // need boolean for current implementation of G4PhysicsVector
    // ==> being phased out
    G4bool isOutRange;
    G4double CAImin = CerenkovAngleIntegrals->GetValue(Pmin, isOutRange);
    ge = CAImax - CAImin;

    if (verboseLevel > 0) {
      G4cout << "CAImin = " << CAImin << G4endl;
      G4cout << "ge = " << ge << G4endl;
    }
  }

  // Calculate number of photons
  G4double NumPhotons = Rfact * charge / eplus * charge / eplus *
                        (dp - ge * BetaInverse * BetaInverse);

  return NumPhotons;
}
