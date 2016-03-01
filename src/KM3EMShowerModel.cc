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
#include "KM3EMShowerModel.hh"

#include "Randomize.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4VProcess.hh"
#include "G4TransportationManager.hh"
#include "KM3TrackInformation.hh"

KM3EMShowerModel::KM3EMShowerModel(std::string modelName, G4Region *envelope)
    : G4VFastSimulationModel(modelName, envelope) {
  EnergyThreshold =
      31.6 * MeV; // this just for speed. Param. slows sim below that threshold
}

KM3EMShowerModel::KM3EMShowerModel(std::string modelName)
    : G4VFastSimulationModel(modelName) {
  EnergyThreshold =
      31.6 * MeV; // this just for speed. Param. slows sim below that threshold
}

KM3EMShowerModel::~KM3EMShowerModel() {
#ifndef G4DISABLE_PARAMETRIZATION
  delete myFlux;
#else
  ;
#endif
}

void KM3EMShowerModel::InitializeFlux(char *infileParam,
                                      double Quantum_Efficiency,
                                      double TotCathodArea) {
#ifndef G4DISABLE_PARAMETRIZATION
  myFlux = new KM3EMEnergyFlux(infileParam, Quantum_Efficiency, TotCathodArea,
                               8, 100.0);
#else
  ;
#endif
}

// the following is for em cascades only.
// do not need hadronic- they are pretty rare and they produce muons that are
// not descibed by the parametrization
// due to the long range.
bool
KM3EMShowerModel::IsApplicable(const G4ParticleDefinition &particleType) {
  if (&particleType == G4Electron::ElectronDefinition() ||
      &particleType == G4Positron::PositronDefinition() ||
      &particleType == G4Gamma::GammaDefinition() // do not to need the pi0 it
                                                  // is discribed by gamma
                                                  // (through decay)
      ) {
    return true;
  }

  else {
    return false;
  }
}

bool KM3EMShowerModel::ModelTrigger(const G4FastTrack &fastTrack) {
  // Applies the parameterisation only if the particle is in Water and the
  // energy is supported by the flux (interpolation model)
  std::string materialName;
  double kineticEnergy;
  materialName = fastTrack.GetPrimaryTrack()->GetMaterial()->GetName();
  kineticEnergy = fastTrack.GetPrimaryTrack()->GetKineticEnergy();
  return ((kineticEnergy > EnergyThreshold) &&
          myFlux->ModelTrigger(kineticEnergy) && (materialName == "Water"));
}

void KM3EMShowerModel::DoIt(const G4FastTrack &fastTrack,
                            G4FastStep &fastStep) {
#if !defined(G4MYEM_PARAMETERIZATION) && !defined(G4MYHA_PARAMETERIZATION)
  static int ooo = 0;
  if (ooo == 0) {
    ooo = 1;
    G4Material *aMaterial = G4Material::GetMaterial("Cathod");
    double MaxQE = -1;
    double PhEneAtMaxQE;
    G4MaterialPropertyVector *aPropertyVector =
        aMaterial->GetMaterialPropertiesTable()->GetProperty("Q_EFF");

    for (size_t i = 0; i < aPropertyVector->GetVectorLength(); i++) {
      double ThisQE = (*aPropertyVector)[i];
      double ThisPhEne = aPropertyVector->Energy(i);
      if (ThisQE > MaxQE) {
        MaxQE = ThisQE;
        PhEneAtMaxQE = ThisPhEne;
      }
    }
    aMaterial = G4Material::GetMaterial("Water");
    G4MaterialPropertyVector *GroupVel =
        aMaterial->GetMaterialPropertiesTable()->GetProperty("GROUPVEL");
    thespeedmaxQE = GroupVel->Value(PhEneAtMaxQE); // coresponds to the maximum
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
  double primaryShowerEnergy =
      fastTrack.GetPrimaryTrack()->GetKineticEnergy();
  // axis of the shower, in global reference frame (normalized to unity):
  G4ThreeVector primaryShowerAxis =
      fastTrack.GetPrimaryTrack()->GetMomentumDirection();
  // starting point of the shower:
  G4ThreeVector primaryShowerPosition =
      fastTrack.GetPrimaryTrack()->GetPosition();
  int idbeam =
      fastTrack.GetPrimaryTrack()->GetParticleDefinition()->GetPDGEncoding();
  // the time of the primary track
  double primaryShowerTime = fastTrack.GetPrimaryTrack()->GetGlobalTime();
  // next we transfer forwards the position a little in order to account for the
  // distance from the particle creation
  // to the maximum of the photon emission (see also generation action in the EM
  // parametrization section)
  /////////////////////////////////////////////////////////////////////////////////////
  // in the next temporary code, we add the Dgamma/c and substract the
  // Delectron/c when we have gamma, since
  // in the parametrization generation we have not shifted the vertex time to
  // -Delectron/c,
  // so the timing information from the tables corresponds to Delectron/c time,
  // if we suppose that
  // all photons are emitted from the shower maximum
  // Normally we should shift the time in parametrization generation by
  // -Delectron/c and here
  // add just Dgamma/c for gamma and Delectron/c for electron,positron.
  double enepos = log10(primaryShowerEnergy / GeV);
  double zpos, zpos11, zpos22;
  if (enepos < 0) {
    zpos11 = 0.90653 + enepos * (0.96676 + enepos * 0.26983);
    zpos22 = 1.22450 + enepos * (0.82081 + enepos * 0.20987);
  } else {
    zpos11 = 0.90132 + enepos * 0.81286;
    zpos22 = 1.20880 + enepos * 0.78600;
  }
  zpos11 *= m;
  zpos22 *= m;
  if (idbeam == 11 || idbeam == -11) {
    zpos = zpos11;
  } else if (idbeam == 22) {
    zpos = zpos22;
    primaryShowerTime += (zpos22 - zpos11) / c_light;
  }
  primaryShowerPosition += zpos * primaryShowerAxis;

  const G4VProcess *theProcess =
      fastTrack.GetPrimaryTrack()->GetCreatorProcess();
  int originalTrackCreatorProcess;
  int originalParentID;
#ifdef G4TRACK_INFORMATION
  if (theProcess != NULL) {
    KM3TrackInformation *info =
        (KM3TrackInformation *)(fastTrack.GetPrimaryTrack()
                                    ->GetUserInformation());
    originalParentID = info->GetOriginalParentID();
    std::string creator = info->GetOriginalTrackCreatorProcess();
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
    originalTrackCreatorProcess = 6;
  }
#else
  originalParentID = 1;
  originalTrackCreatorProcess = 0;
#endif
  int originalInfo =
      (originalParentID - 1) * 10 + originalTrackCreatorProcess;

  /*
  if(theProcess != NULL){
    G4cout << "KM3EMShowerModel::DoIt particle " <<
  fastTrack.GetPrimaryTrack()->GetDefinition()->GetParticleName()
           << " and Energy " <<primaryShowerEnergy<<" MeV "
           <<" created by "<< theProcess->GetProcessName() << G4endl;
  }
  else{
    G4cout << "KM3EMShowerModel::DoIt particle " <<
  fastTrack.GetPrimaryTrack()->GetDefinition()->GetParticleName()
           << " and Energy " <<primaryShowerEnergy<<" MeV "
           <<" which is one of the primary particles"<< G4endl;
  }
  */
  // initialize flux generator for this shower
  static double MaxAbsDist2 =
      myStDetector->MaxAbsDist * myStDetector->MaxAbsDist;
  //  int   PhotonsSurviving=0;
  static int TotalNumberOfTowers = myStDetector->allTowers->size();
  for (int it = 0; it < TotalNumberOfTowers; it++) {
    double dx = (*(myStDetector->allTowers))[it]->position[0] -
                  primaryShowerPosition[0];
    double dy = (*(myStDetector->allTowers))[it]->position[1] -
                  primaryShowerPosition[1];
    double distancetower2 = dx * dx + dy * dy;
    if (distancetower2 < MaxAbsDist2) {
      int TotalNumberOfOMs =
          (*(myStDetector->allTowers))[it]->BenthosIDs->size();
      for (int iot = 0; iot < TotalNumberOfOMs; iot++) {
        int io = (*(*(myStDetector->allTowers))[it]->BenthosIDs)[iot];
        G4ThreeVector FromGeneToOM =
            (*myStDetector->allOMs)[io]->position - primaryShowerPosition;
        double distancein = FromGeneToOM.mag2();
        if (distancein < MaxAbsDist2) {
          distancein = sqrt(distancein);
          FromGeneToOM /= distancein;
          double anglein = primaryShowerAxis.dot(FromGeneToOM);
          myFlux->FindBins(primaryShowerEnergy, distancein, anglein);
          int NumberOfSamples = myFlux->GetNumberOfSamples();
          int icstart, icstop;
          double theFastTime;
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
          for (int isa = 0; isa < NumberOfSamples; isa++) {
            onePE aPE = myFlux->GetSamplePoint();
            //	    G4cout << "OutFromParam "<<distancein/m<<" "<<anglein<<"
            //"<<aPE.costh<<" "<<aPE.phi<<G4endl;  //tempo
            double costh = aPE.costh;
            double sinth = sqrt(1.0 - costh * costh);
            double cosphi = cos(aPE.phi);
            double sinphi = sin(aPE.phi);
            // short	    G4ThreeVector
            // photonDirection=-(sinth*(cosphi*x+sinphi*y)+costh*z);
            G4ThreeVector photonDirection =
                (sinth * (cosphi * x + sinphi * y) + costh * z);
            // short	    double
            // angleThetaDirection=photonDirection.theta();
            // short	    double anglePhiDirection=photonDirection.phi();
            // short	    angleThetaDirection *= 180./M_PI;
            // short	    anglePhiDirection *= 180./M_PI;
            // short	    if(anglePhiDirection < 0.0)anglePhiDirection +=
            // 360.0;
            // short	    int
            // angleDirection=(int)(nearbyint(angleThetaDirection)*1000.0 +
            // nearbyint(anglePhiDirection));
            int ic = int(icstart + (icstop - icstart) * G4UniformRand());
            // short
            // aMySD->InsertExternalHit(ic,theFastTime+aPE.time,originalInfo,angleDirection,-999);
            aMySD->InsertExternalHit(ic, (*myStDetector->allOMs)[io]->position,
                                     theFastTime + aPE.time, originalInfo,
                                     photonDirection);
            //	    PhotonsSurviving++;
          } // for(int isa=0 ; isa<NumberOfSamples ; isa++)
        } // if(distancein<MaxAbsDist2)
      } // for(int io=0;io<TotalNumberOfOMs;io++)
    } // if(distancetower2<MaxAbsDist2)
  } // for(int it=0;it<TotalNumberOfTowers;it++)
//  G4cout << "Total photons created "<< PhotonsSurviving <<G4endl;
#else
  ;
#endif
}