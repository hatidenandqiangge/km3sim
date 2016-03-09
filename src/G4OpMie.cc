#include "G4ios.hh"
#include "G4ExceptionHandler.hh"
#include <stdio.h>
#include "G4OpMie.hh"
#if (defined(G4MYLASER_PARAMETERIZATION) && defined(G4TRACK_INFORMATION)) || \
    (!defined(G4DISABLE_PARAMETRIZATION) &&                                  \
     defined(G4TRACK_INFORMATION))  // newmie
#include "KM3TrackInformation.hh"
#endif

G4OpMie::G4OpMie(const G4String &processName, G4ProcessType type)
    : G4VDiscreteProcess(processName, type) {
  if (verboseLevel > 0) {
    G4cout << GetProcessName() << " is created " << G4endl;
  }
  thePhaseFactors = NULL;
  BuildThePhysicsTable();
}

G4OpMie::~G4OpMie() {
  for (size_t i = 0; i < thePhaseFactors->size(); i++) {
    free((*(thePhaseFactors))[i]);
  }
  thePhaseFactors->clear();
  delete thePhaseFactors;
  delete PhaseRand;
}

G4VParticleChange *G4OpMie::PostStepDoIt(const G4Track &aTrack,
                                         const G4Step &aStep) {
  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle *aParticle = aTrack.GetDynamicParticle();

  if (verboseLevel > 0) {
    G4cout << "Scattering Photon!" << G4endl;
  }

  // find polar angle of new direction w.r.t. old direction

  G4double CosTheta = std::cos(SampleAngle());
  // G4double CosTheta=std::cos(PhaseRand->shoot()*pi);
  G4double SinTheta = std::sqrt(1. - CosTheta * CosTheta);

  // find azimuthal angle of new direction w.r.t. old direction
  G4double rand = G4UniformRand();
  G4double Phi = twopi * rand;
  G4double SinPhi = std::sin(Phi);
  G4double CosPhi = std::cos(Phi);

  G4double unit_x = SinTheta * CosPhi;
  G4double unit_y = SinTheta * SinPhi;
  G4double unit_z = CosTheta;

  G4ThreeVector NewMomentumDirection(unit_x, unit_y, unit_z);
  G4ThreeVector OldMomentumDirection = aParticle->GetMomentumDirection();
  NewMomentumDirection.rotateUz(OldMomentumDirection);

  // find new polarization
  G4ThreeVector OldPolarization = aParticle->GetPolarization();
  G4ThreeVector NewPolarization =
      OldPolarization -
      OldPolarization.dot(NewMomentumDirection) * NewMomentumDirection;
  if (G4UniformRand() < 0.5) NewPolarization = -NewPolarization;
  NewPolarization = NewPolarization.unit();

  aParticleChange.ProposePolarization(NewPolarization);

  aParticleChange.ProposeMomentumDirection(NewMomentumDirection);

  if (verboseLevel > 0) {
    printf(
        "OldNew All %13.6le %13.6le %13.6le %13.6le %13.6le %13.6le %13.6le "
        "%13.6le %13.6le %13.6le %13.6le %13.6le\n",
        OldMomentumDirection.x(), OldMomentumDirection.y(),
        OldMomentumDirection.z(), NewMomentumDirection.x(),
        NewMomentumDirection.y(), NewMomentumDirection.z(), OldPolarization.x(),
        OldPolarization.y(), OldPolarization.z(), NewPolarization.x(),
        NewPolarization.y(), NewPolarization.z());
  }
#if defined(G4MYLASER_PARAMETERIZATION) && defined(G4TRACK_INFORMATION)
  KM3TrackInformation *info =
      (KM3TrackInformation *)(aTrack.GetUserInformation());
  info->KeepScatteringPosition(aTrack.GetPosition(), CosTheta);
#endif
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// BuildThePhysicsTable for the Mie Scattering process (reads only the phase
// function model)
void G4OpMie::BuildThePhysicsTable() {
  // first load Phase factors
  G4double c0, c1, c2, c3, c4, c5, c6;
  if (!thePhaseFactors) {
    thePhaseFactors = new std::vector<PhaseFactors *>;
    FILE *infilePF;
    if ((infilePF = fopen("MiePhaseFactors.in", "r")) == NULL) {
      G4Exception("Error open input Mie Phase Factors file\n", "",
                  FatalException, "");
    } else {
      for (G4int i = 0; i < 56; i++) {  // read the phase functions from
                                        // measurements in shallow waters
        fscanf(infilePF, "%lf %lf %lf %lf %lf %lf %lf\n", &c0, &c1, &c2, &c3,
               &c4, &c5, &c6);
        PhaseFactors *aPhaseFactors =
            (PhaseFactors *)malloc(sizeof(PhaseFactors));
        aPhaseFactors->c0 = 1.0 / c0;
        aPhaseFactors->c1 = c1;
        aPhaseFactors->c2 = c2;
        aPhaseFactors->c3 = c3;
        aPhaseFactors->c4 = c4;
        aPhaseFactors->c5 = c5;
        aPhaseFactors->c6 = c6;
        thePhaseFactors->push_back(aPhaseFactors);
      }
      // read the phase functions parametrized in deep water
      // first one is f4 model. c1 is parameter (a), c0 is the maximum of
      // sin(th)*Ph(th)
      // Ph(th)=1/(1+a^2-2*a*cos(th))^(3/2) unnormalized
      fscanf(infilePF, "%lf %lf\n", &c0, &c1);
      PhaseFactors *aPhaseFactors =
          (PhaseFactors *)malloc(sizeof(PhaseFactors));
      aPhaseFactors->c0 = 1.0 / c0;
      aPhaseFactors->c1 = c1;
      thePhaseFactors->push_back(aPhaseFactors);

      // second one is p0.0075 model. c1 is parameter (p), the rayleigh
      // scattering contribution
      // c2 is parameter (a) of rayleigh contribution, c3 parameter (b) of
      // rayleigh g(th)=a*(1+b*cos^2(th))
      // c4 is parameter (a) of particulate scattering.
      // PhTot(th)=p*g(th)+(1-p)*(1/4pi)*(1-a^2)*Ph(th)
      // c0 is the maximum of sin(th)*PhTot(th)
      fscanf(infilePF, "%lf %lf %lf %lf %lf\n", &c0, &c1, &c2, &c3, &c4);
      aPhaseFactors = (PhaseFactors *)malloc(sizeof(PhaseFactors));
      aPhaseFactors->c0 = 1.0 / c0;
      aPhaseFactors->c1 = c1;
      aPhaseFactors->c2 = c2;
      aPhaseFactors->c3 = c3;
      aPhaseFactors->c4 = c4;
      thePhaseFactors->push_back(aPhaseFactors);

      fclose(infilePF);
    }
  }

  // loop for materials
  const G4MaterialTable *theMaterialTable = G4Material::GetMaterialTable();
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();

  for (G4int i = 0; i < numOfMaterials; i++) {
    if ((*theMaterialTable)[i]->GetName() == "Water") {
      G4MaterialPropertiesTable *aMaterialPropertiesTable =
          (*theMaterialTable)[i]->GetMaterialPropertiesTable();

      if (aMaterialPropertiesTable) {
        G4double DIndexPhaseFunction = 0.0;
        DIndexPhaseFunction =
            aMaterialPropertiesTable->GetConstProperty("MIEPHASE");
        IndexPhaseFunction = (G4int)DIndexPhaseFunction - 1;
      }
    }
  }
  // next build the array, step 0.1 degrees
  for (G4int i = 0; i <= 1800; i++) {
    G4double angle = G4double(i) / 10.0;
    PhaseArray[i] = PhaseFunction(angle * degree);
  }
  PhaseRand = new CLHEP::RandGeneral(PhaseArray, 1801);
}

G4double G4OpMie::SampleAngle(void) {
  G4double angle;
  G4double angleval, randomval;
  do {
    angle = pi * G4UniformRand();
    angleval = PhaseFunction(angle);
    randomval = G4UniformRand();
  } while (angleval < randomval);
  return angle;
}

// phase function*sin(theta) with maximum 1
G4double G4OpMie::PhaseFunction(G4double angle) {
  G4double val;
  G4double x;
  G4double c0, c1, c2, c3, c4, c5, c6;
  if (IndexPhaseFunction < 56) {
    x = sqrt(angle);
    c0 = (*(thePhaseFactors))[IndexPhaseFunction]->c0;
    c1 = (*(thePhaseFactors))[IndexPhaseFunction]->c1;
    c2 = (*(thePhaseFactors))[IndexPhaseFunction]->c2;
    c3 = (*(thePhaseFactors))[IndexPhaseFunction]->c3;
    c4 = (*(thePhaseFactors))[IndexPhaseFunction]->c4;
    c5 = (*(thePhaseFactors))[IndexPhaseFunction]->c5;
    c6 = (*(thePhaseFactors))[IndexPhaseFunction]->c6;
    val = c0 * std::sin(angle) *
          exp(x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * (c5 + x * c6))))));
  } else if (IndexPhaseFunction == 56) {  // f4 model
    c0 = (*(thePhaseFactors))[IndexPhaseFunction]->c0;
    c1 = (*(thePhaseFactors))[IndexPhaseFunction]->c1;
    val = c0 * std::sin(angle) /
          std::pow(1 + c1 * c1 - 2 * c1 * std::cos(angle), 1.5);
  } else if (IndexPhaseFunction == 57) {  // p0.0075 model
    c0 = (*(thePhaseFactors))[IndexPhaseFunction]->c0;
    c1 = (*(thePhaseFactors))[IndexPhaseFunction]->c1;
    c2 = (*(thePhaseFactors))[IndexPhaseFunction]->c2;
    c3 = (*(thePhaseFactors))[IndexPhaseFunction]->c3;
    c4 = (*(thePhaseFactors))[IndexPhaseFunction]->c4;
    val = c0 * std::sin(angle) *
          (c1 * c2 * (1 + c3 * std::pow(cos(angle), 2.0)) +
           (1 - c1) * (1 / (4 * pi)) * (1 - c4 * c4) /
               std::pow(1 + c4 * c4 - 2 * c4 * std::cos(angle), 1.5));
  }
  return val;
}

G4double G4OpMie::GetMeanFreePath(const G4Track &aTrack, G4double,
                                  G4ForceCondition *) {
  const G4DynamicParticle *aParticle = aTrack.GetDynamicParticle();
  const G4Material *aMaterial = aTrack.GetMaterial();

  G4double thePhotonMomentum = aParticle->GetTotalMomentum();

  G4double AttenuationLength = DBL_MAX;

  G4MaterialPropertiesTable *aMaterialPropertyTable =
      aMaterial->GetMaterialPropertiesTable();

  if (aMaterialPropertyTable) {
    G4MaterialPropertyVector *AttenuationLengthVector =
        aMaterialPropertyTable->GetProperty("MIELENGTH");
    if (AttenuationLengthVector) {
      AttenuationLength = AttenuationLengthVector->Value(thePhotonMomentum);
    } else {
      //               G4cout << "No Mie scattering length specified" << G4endl;
    }
  } else {
    //             G4cout << "No Mie scattering length specified" << G4endl;
  }

  return AttenuationLength;
}
