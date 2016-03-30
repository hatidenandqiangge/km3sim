#include "KM3TrackInformation.h"
#include "G4ParticleDefinition.h"
#include "G4VProcess.h"
#include "G4ios.h"

G4Allocator<KM3TrackInformation> aTrackInformationAllocator;

KM3TrackInformation::KM3TrackInformation() {
  originalTrackCreatorProcess = "";
  originalparticleName = "";
  originalEnergy = 0.;
  originalParentID = 0;
  EmittedAsScattered = true;  // newmie
}

KM3TrackInformation::~KM3TrackInformation() { ; }

KM3TrackInformation::KM3TrackInformation(const G4Track *aTrack) {
  originalTrackCreatorProcess = aTrack->GetCreatorProcess()->GetProcessName();
  originalparticleName = aTrack->GetDefinition()->GetParticleName();
  originalEnergy = aTrack->GetTotalEnergy();
  originalParentID = aTrack->GetParentID();
  EmittedAsScattered = false;  // newmie
}
void KM3TrackInformation::SetMoreInformation(const G4Track *aTrack) {
  originalTrackCreatorProcess = aTrack->GetCreatorProcess()->GetProcessName();
  originalparticleName = aTrack->GetDefinition()->GetParticleName();
  originalEnergy = aTrack->GetTotalEnergy();
  originalParentID = aTrack->GetParentID();
}
void KM3TrackInformation::SetMoreInformation(
    const KM3TrackInformation *aTrackInfo) {
  originalTrackCreatorProcess = aTrackInfo->originalTrackCreatorProcess;
  originalparticleName = aTrackInfo->originalparticleName;
  originalEnergy = aTrackInfo->originalEnergy;
  originalParentID = aTrackInfo->originalParentID;
}
KM3TrackInformation::KM3TrackInformation(
    const KM3TrackInformation *aTrackInfo) {
  originalTrackCreatorProcess = aTrackInfo->originalTrackCreatorProcess;
  originalparticleName = aTrackInfo->originalparticleName;
  originalEnergy = aTrackInfo->originalEnergy;
  originalParentID = aTrackInfo->originalParentID;
  EmittedAsScattered = aTrackInfo->EmittedAsScattered;  // newmie
}

void KM3TrackInformation::Print() const {
  G4cout << "Original track Creator Process " << originalTrackCreatorProcess
         << " of particle " << originalparticleName << " with energy "
         << originalEnergy << " and Parent ID " << originalParentID << G4endl;
}
