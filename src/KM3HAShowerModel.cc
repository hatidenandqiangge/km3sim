//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: ExN05EMShowerModel.cc,v 1.9 2004/11/25 23:35:16 mverderi Exp $
// GEANT4 tag $Name: geant4-07-00-patch-01 $
//
#include "KM3HAShowerModel.hh"

#include "Randomize.hh"

#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZeroLong.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Neutron.hh"
#include "G4AntiNeutron.hh"

#include "G4VProcess.hh"
#include "G4TransportationManager.hh"
#include "KM3TrackInformation.hh"

KM3HAShowerModel::KM3HAShowerModel(G4String modelName, G4Region *envelope)
    : G4VFastSimulationModel(modelName, envelope) {
  EnergyMin = 10.0 * GeV;
  EnergyMax = 100.0 * TeV;
}

KM3HAShowerModel::KM3HAShowerModel(G4String modelName)
    : G4VFastSimulationModel(modelName) {
  EnergyMin = 10.0 * GeV;
  EnergyMax = 100.0 * TeV;
}

KM3HAShowerModel::~KM3HAShowerModel() {
#ifndef G4DISABLE_PARAMETRIZATION
  delete myFlux;
#else
  ;
#endif
}

void KM3HAShowerModel::InitializeFlux(char *infileParam,
                                      G4double Quantum_Efficiency,
                                      G4double TotCathodArea) {
#ifndef G4DISABLE_PARAMETRIZATION
  myFlux = new KM3HAEnergyFlux(infileParam, Quantum_Efficiency, TotCathodArea,
                               EnergyMin, EnergyMax);
#else
  ;
#endif
}

// the following is for ha cascades only for pion plus and pion minus generated
// it is applied also too all other frequently generated hadronic particles
G4bool KM3HAShowerModel::IsApplicable(
    const G4ParticleDefinition &particleType) {
  if (&particleType == G4PionPlus::PionPlusDefinition() ||
      &particleType == G4PionMinus::PionMinusDefinition() ||
      &particleType == G4KaonPlus::KaonPlusDefinition() ||
      &particleType == G4KaonMinus::KaonMinusDefinition() ||
      &particleType == G4KaonZeroLong::KaonZeroLongDefinition() ||
      &particleType == G4Proton::ProtonDefinition() ||
      &particleType == G4AntiProton::AntiProtonDefinition() ||
      &particleType == G4Neutron::NeutronDefinition() ||
      &particleType == G4AntiNeutron::AntiNeutronDefinition()) {
    return true;
  }

  else {
    return false;
  }
}

G4bool KM3HAShowerModel::ModelTrigger(const G4FastTrack &fastTrack) {
  // Applies the parameterisation only if the particle is in Water and the
  // energy is supported by the flux (interpolation model)
  G4String materialName;
  G4double Energy;
  materialName = fastTrack.GetPrimaryTrack()->GetMaterial()->GetName();
  //  if(fastTrack.GetPrimaryTrack()->GetDefinition()->GetBaryonNumber() != 0)
  //  //if it is a baryon use kinetic energy
  //    Energy = fastTrack.GetPrimaryTrack()->GetKineticEnergy();
  //  else Energy = fastTrack.GetPrimaryTrack()->GetTotalEnergy(); //if it is a
  //  meson use total energy
  Energy = (fastTrack.GetPrimaryTrack()->GetMomentum())
               .mag();  // in fact the parametrization creation is according to
                        // momentum
  return ((Energy >= EnergyMin) && (Energy <= EnergyMax) &&
          (materialName == "Water"));
}

void KM3HAShowerModel::DoIt(const G4FastTrack &fastTrack,
                            G4FastStep &fastStep) {
#if !defined(G4MYEM_PARAMETERIZATION) && !defined(G4MYHA_PARAMETERIZATION)
  static G4int ooo = 0;
  if (ooo == 0) {
    ooo = 1;
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
    // thespeedmaxQE=GroupVel->GetProperty(3.102*eV); //coresponds to 400nm
  }

  // Kill the parameterised particle:
  fastStep.KillPrimaryTrack();
  fastStep.ProposePrimaryTrackPathLength(0.0);
  fastStep.ProposeTotalEnergyDeposited(
      fastTrack.GetPrimaryTrack()->GetKineticEnergy());
  //  if(fastTrack.GetPrimaryTrack()->GetKineticEnergy() < 10.0*GeV){  //test
  //    fastStep.SetNumberOfSecondaryTracks(0);
  //    return;
  //  }
  // create the optical photons of the shower:
  //-----------------------------------------------------------------------
  // find the parametres of the Cherenkov photons that are going to hit a storey
  //-----------------------------------------------------------------------

  //--kinetic energy of the particle initiating the shower in MeV-----
  G4double primaryShowerEnergy;
  //  if(fastTrack.GetPrimaryTrack()->GetDefinition()->GetBaryonNumber() != 0)
  //  //if it is a baryon use kinetic energy
  //    primaryShowerEnergy = fastTrack.GetPrimaryTrack()->GetKineticEnergy();
  //  else primaryShowerEnergy = fastTrack.GetPrimaryTrack()->GetTotalEnergy();
  //  //if it is a meson use total energy
  primaryShowerEnergy = (fastTrack.GetPrimaryTrack()->GetMomentum())
                            .mag();  // in fact the parametrization creation is
                                     // according to momentum
  // axis of the shower, in global reference frame (normalized to unity):
  G4ThreeVector primaryShowerAxis =
      fastTrack.GetPrimaryTrack()->GetMomentumDirection();
  // starting point of the shower:
  G4ThreeVector primaryShowerPosition =
      fastTrack.GetPrimaryTrack()->GetPosition();
  // next we transfer forwards the position a little in order to account for the
  // distance from the particle creation
  // to the maximum of the photon emission (see also generation action in the EM
  // and HA parametrization section)
  // the shift is fixed for 100GeV pion and kaon zero long
  G4double zpos = 3.46;
  primaryShowerPosition += zpos * m * primaryShowerAxis;
  // the time of the primary track
  G4double primaryShowerTime = fastTrack.GetPrimaryTrack()->GetGlobalTime();
  // the type of particle
  G4int idbeam =
      fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetPDGEncoding();
  // clear stacks

  const G4VProcess *theProcess =
      fastTrack.GetPrimaryTrack()->GetCreatorProcess();
  G4int originalTrackCreatorProcess;
  G4int originalParentID;
#ifdef G4TRACK_INFORMATION
  if (theProcess != NULL) {
    KM3TrackInformation *info =
        (KM3TrackInformation *)(fastTrack.GetPrimaryTrack()
                                    ->GetUserInformation());
    originalParentID = info->GetOriginalParentID();
    G4String creator = info->GetOriginalTrackCreatorProcess();
    if (creator == "KM3Cherenkov")
      originalTrackCreatorProcess = 0;
    else if (creator == "muPairProd")
      originalTrackCreatorProcess = 1;
    else if (creator == "muIoni")
      originalTrackCreatorProcess = 2;
    else if (creator == "muBrems")
      originalTrackCreatorProcess = 3;
    else if (creator == "muonNuclear")
      originalTrackCreatorProcess = 4;
    else if (creator == "Decay")
      originalTrackCreatorProcess = 8;
    else if (creator == "muMinusCaptureAtRest")
      originalTrackCreatorProcess = 9;
    else
      originalTrackCreatorProcess = 5;
  } else {
    originalParentID = fastTrack.GetPrimaryTrack()->GetTrackID();
    originalTrackCreatorProcess = 7;
  }
#else
  originalParentID = 1;
  originalTrackCreatorProcess = 0;
#endif
  G4int originalInfo =
      (originalParentID - 1) * 10 + originalTrackCreatorProcess;

  /*
  if(theProcess != NULL){
    G4cout << "KM3HAShowerModel::DoIt particle " <<
  fastTrack.GetPrimaryTrack()->GetDefinition()->GetParticleName()
           << " and Energy " <<primaryShowerEnergy<<" MeV "
           <<" created by "<< theProcess->GetProcessName() << G4endl;
  }
  else{
    G4cout << "KM3HAShowerModel::DoIt particle " <<
  fastTrack.GetPrimaryTrack()->GetDefinition()->GetParticleName()
           << " and Energy " <<primaryShowerEnergy<<" MeV "
           <<" which is one of the primary particles"<< G4endl;
  }
  */
  // initialize flux generator for this shower
  static G4double MaxAbsDist2 =
      myStDetector->MaxAbsDist * myStDetector->MaxAbsDist;
  G4int PhotonsSurviving = 0;
  static G4int TotalNumberOfTowers = myStDetector->allTowers->size();
  for (int it = 0; it < TotalNumberOfTowers; it++) {
    G4double dx = (*(myStDetector->allTowers))[it]->position[0] -
                  primaryShowerPosition[0];
    G4double dy = (*(myStDetector->allTowers))[it]->position[1] -
                  primaryShowerPosition[1];
    G4double distancetower2 = dx * dx + dy * dy;
    if (distancetower2 < MaxAbsDist2) {
      G4int TotalNumberOfOMs =
          (*(myStDetector->allTowers))[it]->BenthosIDs->size();
      for (G4int iot = 0; iot < TotalNumberOfOMs; iot++) {
        G4int io = (*(*(myStDetector->allTowers))[it]->BenthosIDs)[iot];
        G4ThreeVector FromGeneToOM =
            (*myStDetector->allOMs)[io]->position - primaryShowerPosition;
        G4double distancein = FromGeneToOM.mag2();
        if (distancein < MaxAbsDist2) {
          distancein = sqrt(distancein);
          FromGeneToOM /= distancein;
          G4double anglein = primaryShowerAxis.dot(FromGeneToOM);
          myFlux->FindBins(idbeam, primaryShowerEnergy, distancein, anglein);
          G4int NumberOfSamples = myFlux->GetNumberOfSamples();
          G4int icstart, icstop;
          G4double theFastTime;
          G4ThreeVector x, y, z;
          if (NumberOfSamples > 0) {
            icstart = (*(*myStDetector->allOMs)[io]->CathodsIDs)[0];
            icstop = 1 +
                     (*(*myStDetector->allOMs)[io]->CathodsIDs)
                         [(*myStDetector->allOMs)[io]->CathodsIDs->size() - 1];
            theFastTime = distancein / thespeedmaxQE + primaryShowerTime;
            z = FromGeneToOM;
            y = primaryShowerAxis.cross(z) / sqrt(1.0 - anglein * anglein);
            x = y.cross(z);
          }
          for (G4int isa = 0; isa < NumberOfSamples; isa++) {
            onePE aPE = myFlux->GetSamplePoint();
            //  G4cout << "OutFromParam "<<distancein<<" "<<anglein<<"
            //"<<aPE.costh<<" "<<aPE.phi<<" "<<aPE.time<<G4endl;  //tempo
            G4double costh = aPE.costh;
            G4double sinth = sqrt(1.0 - costh * costh);
            G4double cosphi = cos(aPE.phi);
            G4double sinphi = sin(aPE.phi);
            // short      G4ThreeVector
            // photonDirection=-(sinth*(cosphi*x+sinphi*y)+costh*z);
            G4ThreeVector photonDirection =
                (sinth * (cosphi * x + sinphi * y) + costh * z);
            // short      G4double
            // angleThetaDirection=photonDirection.theta();
            // short      G4double anglePhiDirection=photonDirection.phi();
            // short      angleThetaDirection *= 180./M_PI;
            // short      anglePhiDirection *= 180./M_PI;
            // short      if(anglePhiDirection < 0.0)anglePhiDirection +=
            // 360.0;
            // short      G4int
            // angleDirection=(G4int)(nearbyint(angleThetaDirection)*1000.0 +
            // nearbyint(anglePhiDirection));
            G4int ic = G4int(icstart + (icstop - icstart) * G4UniformRand());
            // short
            // aMySD->InsertExternalHit(ic,theFastTime+aPE.time,originalInfo,angleDirection,-900);
            aMySD->InsertExternalHit(ic, (*myStDetector->allOMs)[io]->position,
                                     theFastTime + aPE.time, originalInfo,
                                     photonDirection);
            PhotonsSurviving++;
          }  // for(G4int isa=0 ; isa<NumberOfSamples ; isa++)
        }    // if(distancein<MaxAbsDist2)
      }      // for(int io=0;io<TotalNumberOfOMs;io++)
    }        // if(distancetower2<MaxAbsDist2)
  }          // for(int it=0;it<TotalNumberOfTowers;it++)
//  G4cout << "Total photons created "<< PhotonsSurviving <<G4endl;
#else
  ;
#endif
}
