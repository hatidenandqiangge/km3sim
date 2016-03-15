#ifndef KM3Hit_h
#define KM3Hit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Transform3D.hh"


class KM3Hit : public G4VHit {
 public:
  KM3Hit();
  ~KM3Hit();
  KM3Hit(const KM3Hit &);
  const KM3Hit &operator=(const KM3Hit &);
  int operator==(const KM3Hit &) const;

  inline void *operator new(size_t);
  inline void operator delete(void *);

 public:
  void SetCathodId(G4int id) { CathodId = id; };
  void SetTime(G4double tt) { time = tt; };
  void SetoriginalInfo(G4int inf) { originalInfo = inf; };
  void SetMany(G4int im) { IMany = im; };
  // short  void SetangleDirection(G4int an) {angleDirection = an;};
  // short  void SetangleIncident(G4int an) {angleIncident = an;};

  G4double GetTime() { return time; };
  G4int GetCathodId() { return CathodId; };
  G4int GetoriginalInfo() { return originalInfo; };
  G4int GetMany() { return IMany; };
  // short  G4int    GetangleDirection() {return angleDirection;};
  // short  G4int    GetangleIncident() {return angleIncident;};

 private:
  G4int CathodId;
  G4double time;
  G4int originalInfo;
  G4int IMany;
  // short  G4int      angleIncident;
  // short  G4int      angleDirection;
};

typedef G4THitsCollection<KM3Hit> KM3HitsCollection;

extern G4Allocator<KM3Hit> KM3HitAllocator;

inline void *KM3Hit::operator new(size_t) {
  void *aHit;
  aHit = (void *)KM3HitAllocator.MallocSingle();
  return aHit;
}

inline void KM3Hit::operator delete(void *aHit) {
  KM3HitAllocator.FreeSingle((KM3Hit *)aHit);
}

#endif
