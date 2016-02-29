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


#include "KM3TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "KM3TrackInformation.hh"

void KM3TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
#ifdef G4TRACK_INFORMATION
#ifdef G4MYLASER_PARAMETERIZATION    
  if(aTrack->GetUserInformation()==0)
#else
  if(aTrack->GetParentID()>0 && aTrack->GetParentID()<=numofInitialParticles)
#endif
  {
    if(aTrack->GetUserInformation()==0){
      KM3TrackInformation* anInfo = new KM3TrackInformation(aTrack);
      G4Track* theTrack = (G4Track*)aTrack;
      theTrack->SetUserInformation(anInfo);
      //write info on evt file about the muon capture or decay secondaries
      if(useANTARESformat){
	G4String theCreatorProcess=aTrack->GetCreatorProcess()->GetProcessName();
	if((theCreatorProcess == "Decay") || (theCreatorProcess == "muMinusCaptureAtRest") ){
	  G4int trackID=aTrack->GetTrackID();
	  G4int parentID=aTrack->GetParentID();
	  G4ThreeVector pos=aTrack->GetPosition();
	  G4ThreeVector ddd=aTrack->GetMomentumDirection();
	  G4double TotalEnergy=aTrack->GetTotalEnergy();
	  G4double time=aTrack->GetGlobalTime();
	  G4int idPDG=aTrack->GetDefinition()->GetPDGEncoding();
	  TheEVTtoWrite->AddMuonDecaySecondaries(trackID,parentID,
						 pos[0]/m,pos[1]/m,pos[2]/m,
						 ddd[0],ddd[1],ddd[2],
						 TotalEnergy,time,idPDG);
	}
      }
    }
    else{
      //      if(aTrack->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition()){
	KM3TrackInformation* info = (KM3TrackInformation*)(aTrack->GetUserInformation());
	if(info->GetEmittedAsScattered() && info->GetOriginalEnergy()==0.0){
	  info->SetMoreInformation(aTrack);
	}
	//      }
    }
  }
#else
  ;
#endif
}


void KM3TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
//   //tempo
//   if(aTrack->GetParentID()==0){
//     G4ThreeVector pDir=aTrack->GetMomentumDirection();
//     G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
//     if(secondaries){
//       size_t nSeco = secondaries->size();
//       if(nSeco>0){
// 	for(size_t i=0;i<nSeco;i++){ 
// 	  G4String aString=(*secondaries)[i]->GetDefinition()->GetParticleName();
// 	  G4String aCrPr=(*secondaries)[i]->GetCreatorProcess()->GetProcessName();
// 	  if(aString != G4String("opticalphoton")){
// 	    G4int idpart;
// 	    G4int idproc;
// 	    if(aString == G4String("gamma"))idpart=0;
// 	    if(aString == G4String("e-"))idpart=1;
// 	    if(aString == G4String("e+"))idpart=2;
// 	    if(aCrPr == G4String("muIoni"))idproc=0;
// 	    if(aCrPr == G4String("muPairProd"))idproc=1;
// 	    if(aCrPr == G4String("muBrems"))idproc=2;
// 	    G4ThreeVector aMomentum=(*secondaries)[i]->GetMomentum();
// 	    G4double momentum=aMomentum.mag();
// 	    if(momentum > 0){
// 	      G4double costheta=pDir.dot(aMomentum)/momentum;
// 	      printf("FromTracking %d %d %le %le\n",idpart,idproc,momentum,costheta);
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
//   //tempo
#ifdef G4TRACK_INFORMATION
  if(aTrack->GetParentID()>0){
    G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
    if(secondaries){
      KM3TrackInformation* info = (KM3TrackInformation*)(aTrack->GetUserInformation());
      size_t nSeco = secondaries->size();
      if(nSeco>0){
	for(size_t i=0;i<nSeco;i++){ 
	  if((*secondaries)[i]->GetUserInformation()==0){
	    KM3TrackInformation* infoNew = new KM3TrackInformation(info);
	    (*secondaries)[i]->SetUserInformation(infoNew);
	  }
	  else{
	    //		    if((*secondaries)[i]->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition()){
	    KM3TrackInformation* infoNew = (KM3TrackInformation*)((*secondaries)[i]->GetUserInformation());
	    if(infoNew->GetEmittedAsScattered() && infoNew->GetOriginalEnergy()==0.0){
	      infoNew->SetMoreInformation(info);
	    }
	    //		    }
	  }
	}
      }
    }
  }
#else
  ;
#endif
}
