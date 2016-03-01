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
  void SetMoreInformation(const G4Track *aTrack); // newmie
  void SetMoreInformation(const KM3TrackInformation *aTrackInfo); // newmie
  ~KM3TrackInformation();

  inline void *operator new(size_t);
  inline void operator delete(void *aTrackInfo);
  inline int operator==(const KM3TrackInformation &right) const {
    return (this == &right);
  }

  void Print() const;

private:
  std::string originalTrackCreatorProcess;
  std::string originalparticleName;
  double originalEnergy;
  int originalParentID;
  bool EmittedAsScattered; // newmie
#ifdef G4MYLASER_PARAMETERIZATION
  std::vector<G4ThreeVector> *ScatteringPositions;
  std::vector<double> *ScatteringAngles;
#endif

public:
  inline std::string GetOriginalTrackCreatorProcess() const {
    return originalTrackCreatorProcess;
  }
  inline std::string GetOriginalParticleName() const {
    return originalparticleName;
  }
  inline double GetOriginalEnergy() const { return originalEnergy; }
  inline int GetOriginalParentID() const { return originalParentID; }
  inline bool GetEmittedAsScattered() const {
    return EmittedAsScattered;
  } // newmie
#ifdef G4MYLASER_PARAMETERIZATION
  inline void KeepScatteringPosition(const G4ThreeVector &ScatteringPosition,
                                     double ScatteringAngle) {
    ScatteringPositions->push_back(ScatteringPosition);
    ScatteringAngles->push_back(ScatteringAngle);
  };
  inline int GetNumberOfScatters(void) {
    return ScatteringPositions->size();
  };
  inline G4ThreeVector GetScatteringPosition(int is) {
    return (*ScatteringPositions)[is];
  };
  inline double GetScatteringAngle(int is) {
    return (*ScatteringAngles)[is];
  };
#endif
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
