#ifndef KM3TrackInformation_h
#define KM3TrackInformation_h 1

#include "globals.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

class KM3TrackInformation : public G4VUserTrackInformation {
 public:
  KM3TrackInformation();
  KM3TrackInformation(const G4Track *aTrack);
  KM3TrackInformation(const KM3TrackInformation *aTrackInfo);
  void SetMoreInformation(const G4Track *aTrack);                  // newmie
  void SetMoreInformation(const KM3TrackInformation *aTrackInfo);  // newmie
  ~KM3TrackInformation();

  inline void *operator new(size_t);
  inline void operator delete(void *aTrackInfo);
  inline int operator==(const KM3TrackInformation &right) const {
    return (this == &right);
  }

  void Print() const;

 private:
  G4String originalTrackCreatorProcess;
  G4String originalparticleName;
  G4double originalEnergy;
  G4int originalParentID;
  G4bool EmittedAsScattered;  // newmie

 public:
  inline G4String GetOriginalTrackCreatorProcess() const {
    return originalTrackCreatorProcess;
  }
  inline G4String GetOriginalParticleName() const {
    return originalparticleName;
  }
  inline G4double GetOriginalEnergy() const { return originalEnergy; }
  inline G4int GetOriginalParentID() const { return originalParentID; }
  inline G4bool GetEmittedAsScattered() const {
    return EmittedAsScattered;
  }  // newmie
};

extern G4Allocator<KM3TrackInformation> aTrackInformationAllocator;

inline void *KM3TrackInformation::operator new(size_t) {
  void *aTrackInfo;
  aTrackInfo = (void *)aTrackInformationAllocator.MallocSingle();
  return aTrackInfo;
}

inline void KM3TrackInformation::operator delete(void *aTrackInfo) {
  aTrackInformationAllocator.FreeSingle((KM3TrackInformation *)aTrackInfo);
}

#endif
