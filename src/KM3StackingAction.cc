#include "KM3StackingAction.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4VProcess.hh"
#include <math.h>
#include "G4StackManager.hh"
// the following was added to see what initial hadrons can give muons
//#include "KM3TrackInformation.hh"
//#ifdef G4MYHAMUONS_PARAMETERIZATION
#include "G4RunManager.hh"
#include "KM3PrimaryGeneratorAction.hh"
//#endif
#ifdef G4ENABLE_MIE
#ifndef G4DISABLE_PARAMETRIZATION

#include "KM3SD.hh"
#include "G4SDManager.hh"
#include "KM3Cherenkov.hh"
#include "G4ProcessTable.hh"

#endif
#endif

KM3StackingAction::KM3StackingAction() {
  ;
}

KM3StackingAction::~KM3StackingAction() {
  ;
}

//#############################################################################
G4ClassificationOfNewTrack KM3StackingAction::ClassifyNewTrack(
    const G4Track *aTrack) {
  static KM3PrimaryGeneratorAction *aGeneAction =
      (KM3PrimaryGeneratorAction *)(G4RunManager::GetRunManager()
                                        ->GetUserPrimaryGeneratorAction());
  G4double kineticEnergy;
  G4ThreeVector x0;
  G4ThreeVector p0;
  G4ThreeVector distanceV;
  G4double direction, distanceRho2;
  static G4double detectorMaxRho2 =
      MyStDetector->detectorMaxRho * MyStDetector->detectorMaxRho;

  // here kill tracks that have already killed by other classes
  if ((aTrack->GetTrackStatus() == fStopAndKill) ||
      (aTrack->GetTrackStatus() == fKillTrackAndSecondaries))
    return fKill;


#ifdef G4HADRONIC_COMPILE
#ifndef G4DISABLE_PARAMETRIZATION
  // here we kill any generated muon with kinetic energy >1GeV that is not
  // primary particle since we have added this to the genarator action
  if (aTrack->GetDefinition() == G4MuonPlus::MuonPlusDefinition() ||
      aTrack->GetDefinition() == G4MuonMinus::MuonMinusDefinition()) {
    if ((aTrack->GetParentID() != 0) &&
        (aTrack->GetKineticEnergy() > 1.0 * GeV))
      return fKill;
  }
#endif
#endif
  //----------------------------------------------------------------------------------------

  if (aTrack->GetDefinition() !=
      G4OpticalPhoton::OpticalPhotonDefinition()) {  // first check that is not
                                                     // a
                                                     // photon to save time
    kineticEnergy = aTrack->GetKineticEnergy();

    // threshold for electron cerenkov production (not applicable for positron
    // due to anihhilation
    if ((aTrack->GetDefinition()->GetParticleName() == "e-") &&
        (kineticEnergy < 240 * keV))
      return fKill;
    // threshold for gamma potentially absorbed by an atomic electron
    if ((aTrack->GetDefinition()->GetParticleName() == "gamma") &&
        (kineticEnergy < 240 * keV))
      return fKill;
    // kill produced neutrinos
    if ((aTrack->GetDefinition()->GetParticleType() == "lepton") &&
        (aTrack->GetDefinition()->GetPDGMass() == 0.0))
      return fKill;

    x0 = aTrack->GetPosition();
    distanceV = x0 - MyStDetector->detectorCenter;
    distanceRho2 = distanceV[0] * distanceV[0] + distanceV[1] * distanceV[1];

    // if the particle is not muon and is created outside the can kill it
    if (aTrack->GetDefinition() != G4MuonPlus::MuonPlusDefinition() &&
        aTrack->GetDefinition() != G4MuonMinus::MuonMinusDefinition()) {
      if ((x0[2] < MyStDetector->bottomPosition) ||
          (distanceRho2 > detectorMaxRho2) ||
          (x0[2] > MyStDetector->detectorMaxz))
        return fKill;
#ifdef G4ENABLE_MIE
#ifndef G4DISABLE_PARAMETRIZATION
      // delta rays parametrization section. Works only for muons primaries (not
      // for showers ve)
      // kill e- from muIoni and keep position and deposited energy information
      if (aTrack->GetCreatorProcess() !=
          NULL) {  // it is not a primary particle
        G4int parentID = aTrack->GetParentID();
        if ((aTrack->GetCreatorProcess()->GetProcessName() == "muIoni") &&
            parentID <= aGeneAction->numberofParticles &&
            (kineticEnergy < 31.6)) {  // delta rays up to 31.6MeV kinene
          // save point and energy deposited above threshold
          idprikeep->push_back(parentID);
          depenekeep->push_back(kineticEnergy - 0.24);
          poskeep->push_back(x0);
          timekeep->push_back(aTrack->GetGlobalTime());
          //// here we create the photons from direct emission and delta rays
          /// before the end of the event
          //// in order for the memory usage not to grow too large (specially
          /// dealing with multi muon signal - showers)
          //// We also create the these photons at the end of the event through
          /// NewStage (double code)
          if (idprikeep->size() == 500000) CreateAllWaitingPhotons();
          ////
          // before  if(indexkeep == 0){
          // before    return fWaiting; //i do this in order to add the hits
          // of
          // the delta rays at the end of the event
          // before    indexkeep=1;
          // before  }
          // before  else return fKill;
          return fKill;  // now
        }
      }
////////////////////////////////////
#endif
#endif
      return fUrgent;
    } else {  // if it is a muon kill it only if is not going to cross the can

      p0 = aTrack->GetMomentumDirection();
      if ((x0[2] < MyStDetector->bottomPosition) && (p0[2] < 0))
        return fKill;  // goes down while below the can
      if ((x0[2] > MyStDetector->detectorMaxz) && (p0[2] > 0))
        return fKill;  // goes up while above the can
      direction = p0[0] * distanceV[0] + p0[1] * distanceV[1];
      if ((distanceRho2 > detectorMaxRho2) && (direction > 0))
        return fKill;  // goes away while outside the can

      // first check if it is inside the can
      G4double rxy2 = x0[0] * x0[0] + x0[1] * x0[1];
      if ((rxy2 < detectorMaxRho2) && (x0[2] > MyStDetector->bottomPosition) &&
          (x0[2] < MyStDetector->detectorMaxz))
        return fUrgent;
      G4double Ttop = (MyStDetector->detectorMaxz - x0[2]) / p0[2];
      if (Ttop > 0) {
        G4double Xtop = x0[0] + Ttop * p0[0] - MyStDetector->detectorCenter[0];
        G4double Ytop = x0[1] + Ttop * p0[1] - MyStDetector->detectorCenter[1];
        G4double dRhoTop = Xtop * Xtop + Ytop * Ytop;
        if (dRhoTop < detectorMaxRho2) return fUrgent;
      }
      G4double Tbottom = (MyStDetector->bottomPosition - x0[2]) / p0[2];
      if (Tbottom > 0) {
        G4double Xbottom =
            x0[0] + Tbottom * p0[0] - MyStDetector->detectorCenter[0];
        G4double Ybottom =
            x0[1] + Tbottom * p0[1] - MyStDetector->detectorCenter[1];
        G4double dRhoBottom = Xbottom * Xbottom + Ybottom * Ybottom;
        if (dRhoBottom < detectorMaxRho2) return fUrgent;
      }
      G4double a = p0[0] * p0[0] + p0[1] * p0[1];
      G4double b = x0[0] * p0[0] + x0[1] * p0[1];
      G4double c = x0[0] * x0[0] + x0[1] * x0[1] - detectorMaxRho2;
      G4double dia = b * b - a * c;
      if (dia > 0) {
        dia = sqrt(dia);
        G4double SideDist1 = (-b - dia) / a;
        G4double sidez1 = x0[2] + SideDist1 * p0[2];
        if ((sidez1 > MyStDetector->bottomPosition) &&
            (sidez1 < MyStDetector->detectorMaxz) && (SideDist1 > 0))
          return fUrgent;
        G4double SideDist2 = (-b + dia) / a;
        G4double sidez2 = x0[2] + SideDist2 * p0[2];
        if ((sidez2 > MyStDetector->bottomPosition) &&
            (sidez2 < MyStDetector->detectorMaxz) && (SideDist2 > 0))
          return fUrgent;
      }
      return fKill;
    }
  }
  return fUrgent;
}


void KM3StackingAction::NewStage() {
#ifdef G4ENABLE_MIE
#ifndef G4DISABLE_PARAMETRIZATION
  // delta rays parametrization section. Works only for muons primaries (not for
  // showers ve)
  // initialization for the group velocity at maximum qe and MySD pointer
  static G4double thespeedmaxQE = -1.0;
  static KM3SD *aMySD = NULL;
  if (thespeedmaxQE < 0) {
    G4Material *aMaterial = G4Material::GetMaterial("Cathod");
    G4double MaxQE = -1;
    G4double PhEneAtMaxQE;
    G4MaterialPropertyVector *aPropertyVector =
        aMaterial->GetMaterialPropertiesTable()->GetProperty("Q_EFF");
    for (size_t i = 0; i < aPropertyVector->GetVectorLength(); i++) {
      G4double ThisQE = (*aPropertyVector)[i];
      G4double ThisPhEne = aPropertyVector->Energy(i);
      if (ThisQE > MaxQE) {
        MaxQE = ThisQE;
        PhEneAtMaxQE = ThisPhEne;
      }
    }
    aMaterial = G4Material::GetMaterial("Water");
    G4MaterialPropertyVector *GroupVel =
        aMaterial->GetMaterialPropertiesTable()->GetProperty("GROUPVEL");
    thespeedmaxQE = GroupVel->Value(PhEneAtMaxQE);  // coresponds to the maximum
                                                    // qe each time. This is the
                                                    // right one
    G4SDManager *SDman = G4SDManager::GetSDMpointer();
    aMySD = (KM3SD *)SDman->FindSensitiveDetector(G4String("mydetector1/MySD"),
                                                  true);
    if (myFlux == NULL)
      myFlux = new KM3EMDeltaFlux(MyStDetector->EMParametrization_FILE,
                                  MyStDetector->Quantum_Efficiency,
                                  MyStDetector->TotCathodArea);
  }
  // end of initialization
  static G4int originalTrackCreatorProcess = 2;  // it is always muIoni
  static G4int TotalNumberOfTowers = MyStDetector->allTowers->size();
  static G4double MaxAbsDist2 =
      MyStDetector->MaxAbsDist * MyStDetector->MaxAbsDist;
  static G4double distbin2 =
      (100.0 * cm) *
      (100.0 * cm);  // is the distance binning in the delta rays for gathering
  size_t arraysize = idprikeep->size();
  //  G4cout << "----------- size from stacking NewStage -----------
  //  "<<arraysize<<G4endl;
  size_t counter = 0;
  while (counter < arraysize) {
    G4int idpri = (*idprikeep)[counter];
    while ((counter < arraysize) && (idpri == (*idprikeep)[counter])) {
      G4int originalInfo = ((*idprikeep)[counter] - 1) * 10 +
                           originalTrackCreatorProcess;  // to write in MySD
      G4ThreeVector p0(0.0, 0.0, 0.0);
      G4ThreeVector pospri = (*poskeep)[counter];
      G4double depene = 0.0;
      G4ThreeVector thispos(0.0, 0.0, 0.0);
      G4double thistime = 0.0;
      while ((counter < arraysize) && (idpri == (*idprikeep)[counter]) &&
             ((pospri - (*poskeep)[counter]).mag2() < distbin2)) {
        depene += (*depenekeep)[counter];
        thispos += (*depenekeep)[counter] * (*poskeep)[counter];
        thistime += (*depenekeep)[counter] * (*timekeep)[counter];
        counter++;
      }
      G4ThreeVector p0this = (*poskeep)[counter - 1] - pospri;
      if (p0this.mag2() > 0.0) p0 = p0this;
      G4double step = p0.mag();
      if (step > 0.0) {  // discard delta rays that are alone in distance
                         // inetrval (too low signal)
        p0 /= step;
        thispos /= depene;
        thistime /= depene;
        for (int it = 0; it < TotalNumberOfTowers; it++) {
          G4double dx =
              (*(MyStDetector->allTowers))[it]->position[0] - thispos[0];
          G4double dy =
              (*(MyStDetector->allTowers))[it]->position[1] - thispos[1];
          G4double distancetower2 = dx * dx + dy * dy;
          if (distancetower2 < MaxAbsDist2) {
            G4int TotalNumberOfOMs =
                (*(MyStDetector->allTowers))[it]->BenthosIDs->size();
            for (int iot = 0; iot < TotalNumberOfOMs; iot++) {
              G4int io = (*(*(MyStDetector->allTowers))[it]->BenthosIDs)[iot];
              G4ThreeVector FromGeneToOM =
                  (*MyStDetector->allOMs)[io]->position - thispos;
              G4double distancein = FromGeneToOM.mag2();
              if (distancein < MaxAbsDist2) {
                distancein = sqrt(distancein);
                FromGeneToOM /= distancein;
                G4double anglein = p0.dot(FromGeneToOM);
                myFlux->FindBins(depene, distancein, anglein);  // here change
                G4int NumberOfSamples =
                    myFlux->GetNumberOfSamples();  // here change
                G4int icstart, icstop;
                G4double theFastTime;
                G4ThreeVector x, y, z;
                if (NumberOfSamples > 0) {
                  icstart = (*(*MyStDetector->allOMs)[io]->CathodsIDs)[0];
                  icstop =
                      1 +
                      (*(*MyStDetector->allOMs)[io]->CathodsIDs)
                          [(*MyStDetector->allOMs)[io]->CathodsIDs->size() - 1];
                  theFastTime = distancein / thespeedmaxQE + thistime;
                  z = FromGeneToOM;
                  y = p0.cross(z) / sqrt(1.0 - anglein * anglein);
                  x = y.cross(z);
                }
                for (G4int isa = 0; isa < NumberOfSamples; isa++) {
                  onePE aPE = myFlux->GetSamplePoint();
                  G4double costh = aPE.costh;  // here change
                  G4double sinth = sqrt(1.0 - costh * costh);
                  G4double cosphi = cos(aPE.phi);  // here change
                  G4double sinphi = sin(aPE.phi);  // here change
                  // short      G4ThreeVector
                  // photonDirection=-(sinth*(cosphi*x+sinphi*y)+costh*z);
                  G4ThreeVector photonDirection =
                      (sinth * (cosphi * x + sinphi * y) + costh * z);
                  // short      G4double
                  // angleThetaDirection=photonDirection.theta();
                  // short      G4double
                  // anglePhiDirection=photonDirection.phi();
                  // short      angleThetaDirection *= 180./M_PI;
                  // short      anglePhiDirection *= 180./M_PI;
                  // short      if(anglePhiDirection <
                  // 0.0)anglePhiDirection
                  // += 360.0;
                  // short      G4int
                  // angleDirection=(G4int)(nearbyint(angleThetaDirection)*1000.0
                  // + nearbyint(anglePhiDirection));
                  G4int ic =
                      G4int(icstart + (icstop - icstart) * G4UniformRand());
                  // short
                  // aMySD->InsertExternalHit(ic,theFastTime+aPE.time,originalInfo,angleDirection,-997);
                  aMySD->InsertExternalHit(
                      ic, (*MyStDetector->allOMs)[io]->position,
                      theFastTime + aPE.time, originalInfo, photonDirection);
                }  // for(G4int isa=0 ; isa<NumberOfSamples ; isa++){
              }    // if(distancein<MyStDetector->MaxAbsDist){
            }      // for(int io=0;io<TotalNumberOfOMs;io++){
          }        // if(distancetower2<MaxAbsDist2)
        }          // for(int it=0;it<TotalNumberOfTowers;it++)
      }            // if(step > 0.0){
    }              // while( idpri == (*idprikeep)[counter] ){
  }                // while (counter<arraysize){
  stackManager
      ->clear();  // delete all waiting tracks from delta rays and end the event
  ////////////////////////////////////
  // section direct ch-photons
  static KM3Cherenkov *myCher = NULL;
  if (myCher == NULL) {
    myCher = (KM3Cherenkov *)(G4ProcessTable::GetProcessTable()->FindProcess(
        "KM3Cherenkov", "mu-"));
  }
  myCher->CreateDirectPhotons();
///////////////////////////////////
endif
#endif

}

#ifdef G4ENABLE_MIE
#ifndef G4DISABLE_PARAMETRIZATION
void KM3StackingAction::CreateAllWaitingPhotons() {
  // delta rays parametrization section. Works only for muons primaries (not for
  // showers ve)
  // initialization for the group velocity at maximum qe and MySD pointer
  static G4double thespeedmaxQE = -1.0;
  static KM3SD *aMySD = NULL;
  if (thespeedmaxQE < 0) {
    G4Material *aMaterial = G4Material::GetMaterial("Cathod");
    G4double MaxQE = -1;
    G4double PhEneAtMaxQE;
    G4MaterialPropertyVector *aPropertyVector =
        aMaterial->GetMaterialPropertiesTable()->GetProperty("Q_EFF");
    for (size_t i = 0; i < aPropertyVector->GetVectorLength(); i++) {
      G4double ThisQE = (*aPropertyVector)[i];
      G4double ThisPhEne = aPropertyVector->Energy(i);
      if (ThisQE > MaxQE) {
        MaxQE = ThisQE;
        PhEneAtMaxQE = ThisPhEne;
      }
    }
    aMaterial = G4Material::GetMaterial("Water");
    G4MaterialPropertyVector *GroupVel =
        aMaterial->GetMaterialPropertiesTable()->GetProperty("GROUPVEL");
    thespeedmaxQE = GroupVel->Value(PhEneAtMaxQE);  // coresponds to the maximum
                                                    // qe each time. This is the
                                                    // right one
    G4SDManager *SDman = G4SDManager::GetSDMpointer();
    aMySD = (KM3SD *)SDman->FindSensitiveDetector(G4String("mydetector1/MySD"),
                                                  true);
    if (myFlux == NULL)
      myFlux = new KM3EMDeltaFlux(MyStDetector->EMParametrization_FILE,
                                  MyStDetector->Quantum_Efficiency,
                                  MyStDetector->TotCathodArea);
  }
  // end of initialization
  static G4int originalTrackCreatorProcess = 2;  // it is always muIoni
  static G4int TotalNumberOfTowers = MyStDetector->allTowers->size();
  static G4double MaxAbsDist2 =
      MyStDetector->MaxAbsDist * MyStDetector->MaxAbsDist;
  static G4double distbin2 =
      (100.0 * cm) *
      (100.0 * cm);  // is the distance binning in the delta rays for gathering
  size_t arraysize = idprikeep->size();
  //  G4cout << "----------- size from stacking CreateAllWaitingPhotons
  //  ----------- "<<arraysize<<G4endl;
  size_t counter = 0;
  while (counter < arraysize) {
    G4int idpri = (*idprikeep)[counter];
    while ((counter < arraysize) && (idpri == (*idprikeep)[counter])) {
      G4int originalInfo = ((*idprikeep)[counter] - 1) * 10 +
                           originalTrackCreatorProcess;  // to write in MySD
      G4ThreeVector p0(0.0, 0.0, 0.0);
      G4ThreeVector pospri = (*poskeep)[counter];
      G4double depene = 0.0;
      G4ThreeVector thispos(0.0, 0.0, 0.0);
      G4double thistime = 0.0;
      while ((counter < arraysize) && (idpri == (*idprikeep)[counter]) &&
             ((pospri - (*poskeep)[counter]).mag2() < distbin2)) {
        depene += (*depenekeep)[counter];
        thispos += (*depenekeep)[counter] * (*poskeep)[counter];
        thistime += (*depenekeep)[counter] * (*timekeep)[counter];
        counter++;
      }
      G4ThreeVector p0this = (*poskeep)[counter - 1] - pospri;
      if (p0this.mag2() > 0.0) p0 = p0this;
      G4double step = p0.mag();
      if (step >
          0.0) {  // discard delta rays that are alone in distance inerval
                  // (too low signal)
        p0 /= step;
        thispos /= depene;
        thistime /= depene;
        for (int it = 0; it < TotalNumberOfTowers; it++) {
          G4double dx =
              (*(MyStDetector->allTowers))[it]->position[0] - thispos[0];
          G4double dy =
              (*(MyStDetector->allTowers))[it]->position[1] - thispos[1];
          G4double distancetower2 = dx * dx + dy * dy;
          if (distancetower2 < MaxAbsDist2) {
            G4int TotalNumberOfOMs =
                (*(MyStDetector->allTowers))[it]->BenthosIDs->size();
            for (int iot = 0; iot < TotalNumberOfOMs; iot++) {
              G4int io = (*(*(MyStDetector->allTowers))[it]->BenthosIDs)[iot];
              G4ThreeVector FromGeneToOM =
                  (*MyStDetector->allOMs)[io]->position - thispos;
              G4double distancein = FromGeneToOM.mag2();
              if (distancein < MaxAbsDist2) {
                distancein = sqrt(distancein);
                FromGeneToOM /= distancein;
                G4double anglein = p0.dot(FromGeneToOM);
                myFlux->FindBins(depene, distancein, anglein);  // here change
                G4int NumberOfSamples =
                    myFlux->GetNumberOfSamples();  // here change
                G4int icstart, icstop;
                G4double theFastTime;
                G4ThreeVector x, y, z;
                if (NumberOfSamples > 0) {
                  icstart = (*(*MyStDetector->allOMs)[io]->CathodsIDs)[0];
                  icstop =
                      1 +
                      (*(*MyStDetector->allOMs)[io]->CathodsIDs)
                          [(*MyStDetector->allOMs)[io]->CathodsIDs->size() - 1];
                  theFastTime = distancein / thespeedmaxQE + thistime;
                  z = FromGeneToOM;
                  y = p0.cross(z) / sqrt(1.0 - anglein * anglein);
                  x = y.cross(z);
                }
                for (G4int isa = 0; isa < NumberOfSamples; isa++) {
                  onePE aPE = myFlux->GetSamplePoint();
                  G4double costh = aPE.costh;  // here change
                  G4double sinth = sqrt(1.0 - costh * costh);
                  G4double cosphi = cos(aPE.phi);  // here change
                  G4double sinphi = sin(aPE.phi);  // here change
                  // short      G4ThreeVector
                  // photonDirection=-(sinth*(cosphi*x+sinphi*y)+costh*z);
                  G4ThreeVector photonDirection =
                      (sinth * (cosphi * x + sinphi * y) + costh * z);
                  // short      G4double
                  // angleThetaDirection=photonDirection.theta();
                  // short      G4double
                  // anglePhiDirection=photonDirection.phi();
                  // short      angleThetaDirection *= 180./M_PI;
                  // short      anglePhiDirection *= 180./M_PI;
                  // short      if(anglePhiDirection <
                  // 0.0)anglePhiDirection
                  // += 360.0;
                  // short      G4int
                  // angleDirection=(G4int)(nearbyint(angleThetaDirection)*1000.0
                  // + nearbyint(anglePhiDirection));
                  G4int ic =
                      G4int(icstart + (icstop - icstart) * G4UniformRand());
                  // short
                  // aMySD->InsertExternalHit(ic,theFastTime+aPE.time,originalInfo,angleDirection,-997);
                  aMySD->InsertExternalHit(
                      ic, (*MyStDetector->allOMs)[io]->position,
                      theFastTime + aPE.time, originalInfo, photonDirection);
                }  // for(G4int isa=0 ; isa<NumberOfSamples ; isa++){
              }    // if(distancein<MyStDetector->MaxAbsDist){
            }      // for(int io=0;io<TotalNumberOfOMs;io++){
          }        // if(distancetower2<MaxAbsDist2)
        }          // for(int it=0;it<TotalNumberOfTowers;it++)
      }            // if(step > 0.0){
    }              // while( idpri == (*idprikeep)[counter] ){
  }                // while (counter<arraysize){
  idprikeep->clear();
  depenekeep->clear();
  poskeep->clear();
  timekeep->clear();
  ////////////////////////////////////
  // section direct ch-photons
  static KM3Cherenkov *myCher = NULL;
  if (myCher == NULL) {
    myCher = (KM3Cherenkov *)(G4ProcessTable::GetProcessTable()->FindProcess(
        "KM3Cherenkov", "mu-"));
  }
  myCher->CreateDirectPhotons();
  ///////////////////////////////////
}
#endif
#endif
#endif

void KM3StackingAction::PrepareNewEvent() {
#if !defined(G4MYEM_PARAMETERIZATION) && !defined(G4MYHA_PARAMETERIZATION)
#ifdef G4ENABLE_MIE
#ifndef G4DISABLE_PARAMETRIZATION
  // delta rays parametrization section. Works only for muons primaries (not for
  // showers ve)
  // before  indexkeep=0;
  idprikeep->clear();
  depenekeep->clear();
  poskeep->clear();
  timekeep->clear();
////////////////////////////////////
#endif
#endif

}

void KM3StackingAction::SetDetector(KM3Detector *adet) { MyStDetector = adet; }
