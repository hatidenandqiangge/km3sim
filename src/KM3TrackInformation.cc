#include "KM3TrackInformation.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"

G4Allocator<KM3TrackInformation> aTrackInformationAllocator;

KM3TrackInformation::KM3TrackInformation()
{
  originalTrackCreatorProcess = "";
  originalparticleName = "";
  originalEnergy = 0.;
  originalParentID = 0;
  EmittedAsScattered = true; //newmie
}

KM3TrackInformation::~KM3TrackInformation()
{
#ifdef G4MYLASER_PARAMETERIZATION
  ScatteringPositions->clear();
  ScatteringAngles->clear();
  delete ScatteringPositions;
  delete ScatteringAngles;
#endif
  ;
}

KM3TrackInformation::KM3TrackInformation(const G4Track* aTrack)
{
#ifndef G4MYLASER_PARAMETERIZATION
  originalTrackCreatorProcess = aTrack->GetCreatorProcess()->GetProcessName();
#endif
  originalparticleName = aTrack->GetDefinition()->GetParticleName();
  originalEnergy = aTrack->GetTotalEnergy();
  originalParentID = aTrack->GetParentID();
  EmittedAsScattered = false; //newmie
#ifdef G4MYLASER_PARAMETERIZATION
  ScatteringPositions = new std::vector<G4ThreeVector>;
  ScatteringAngles = new std::vector<double>;
  KeepScatteringPosition(aTrack->GetPosition(),1.0);
#endif
}
void KM3TrackInformation::SetMoreInformation(const G4Track* aTrack)
{
  originalTrackCreatorProcess = aTrack->GetCreatorProcess()->GetProcessName();
  originalparticleName = aTrack->GetDefinition()->GetParticleName();
  originalEnergy = aTrack->GetTotalEnergy();
  originalParentID = aTrack->GetParentID();
}
void KM3TrackInformation::SetMoreInformation(const KM3TrackInformation* aTrackInfo)
{
  originalTrackCreatorProcess = aTrackInfo->originalTrackCreatorProcess;
  originalparticleName = aTrackInfo->originalparticleName;
  originalEnergy = aTrackInfo->originalEnergy;
  originalParentID = aTrackInfo->originalParentID;
}
KM3TrackInformation::KM3TrackInformation(const KM3TrackInformation* aTrackInfo)
{
  originalTrackCreatorProcess = aTrackInfo->originalTrackCreatorProcess;
  originalparticleName = aTrackInfo->originalparticleName;
  originalEnergy = aTrackInfo->originalEnergy;
  originalParentID = aTrackInfo->originalParentID;
  EmittedAsScattered = aTrackInfo->EmittedAsScattered; //newmie
}

void KM3TrackInformation::Print() const
{
    G4cout 
     << "Original track Creator Process " << originalTrackCreatorProcess 
     << " of particle " << originalparticleName 
     << " with energy " << originalEnergy  
     << " and Parent ID " << originalParentID << G4endl;
}
