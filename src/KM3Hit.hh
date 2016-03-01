
#ifndef KM3Hit_h
#define KM3Hit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Transform3D.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
  void SetCathodId(int id) { CathodId = id; };
  void SetTime(double tt) { time = tt; };
  void SetoriginalInfo(int inf) { originalInfo = inf; };
  void SetMany(int im) { IMany = im; };
#ifdef G4MYLASER_PARAMETERIZATION
  void SetIManyScatters(int im) { IManyScatters = im; };
  void SetScatteringSteps(int is, double step) {
    ScatteringSteps[is] = step;
  };
  void SetScatteringAngles(int is, double angle) {
    ScatteringAngles[is] = angle;
  };
#endif
  // short  void SetangleDirection(int an) {angleDirection = an;};
  // short  void SetangleIncident(int an) {angleIncident = an;};

  double GetTime() { return time; };
  int GetCathodId() { return CathodId; };
  int GetoriginalInfo() { return originalInfo; };
  int GetMany() { return IMany; };
#ifdef G4MYLASER_PARAMETERIZATION
  int GetIManyScatters(void) { return IManyScatters; };
  double GetScatteringSteps(int is) { return ScatteringSteps[is]; };
  double GetScatteringAngles(int is) { return ScatteringAngles[is]; };
#endif
  // short  int    GetangleDirection() {return angleDirection;};
  // short  int    GetangleIncident() {return angleIncident;};

private:
  int CathodId;
  double time;
  int originalInfo;
  int IMany;
#ifdef G4MYLASER_PARAMETERIZATION
  int IManyScatters;
  double ScatteringSteps[100];
  double ScatteringAngles[100];
#endif
  // short  int      angleIncident;
  // short  int      angleDirection;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<KM3Hit> KM3HitsCollection;

extern G4Allocator<KM3Hit> KM3HitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void *KM3Hit::operator new(size_t) {
  void *aHit;
  aHit = (void *)KM3HitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void KM3Hit::operator delete(void *aHit) {
  KM3HitAllocator.FreeSingle((KM3Hit *)aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
