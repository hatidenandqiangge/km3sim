
#ifndef KM3Hit_h
#define KM3Hit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Transform3D.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class KM3Hit : public G4VHit
{
  public:

  KM3Hit();
  ~KM3Hit();
  KM3Hit(const KM3Hit&);
  const KM3Hit& operator=(const KM3Hit&);
  int operator==(const KM3Hit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  

  public:
  
  void SetCathodId(G4int id) {CathodId=id;};        
  void SetTime      (G4double tt){ time = tt; };
  void SetoriginalInfo(G4int inf) {originalInfo = inf;};
  void SetMany (G4int im) {IMany = im;};
#ifdef G4MYLASER_PARAMETERIZATION
  void SetIManyScatters(G4int im){IManyScatters=im;};
  void SetScatteringSteps(G4int is, G4double step){ScatteringSteps[is]=step;};
  void SetScatteringAngles(G4int is, G4double angle){ScatteringAngles[is]=angle;};
#endif
  //short  void SetangleDirection(G4int an) {angleDirection = an;};
  //short  void SetangleIncident(G4int an) {angleIncident = an;};
      
  G4double GetTime(){ return time; };
  G4int    GetCathodId() {return CathodId;};
  G4int    GetoriginalInfo() {return originalInfo;};
  G4int    GetMany() {return IMany;};
#ifdef G4MYLASER_PARAMETERIZATION
  G4int      GetIManyScatters(void){return IManyScatters;};
  G4double   GetScatteringSteps(G4int is){return ScatteringSteps[is];};
  G4double   GetScatteringAngles(G4int is){ return ScatteringAngles[is];};
#endif
  //short  G4int    GetangleDirection() {return angleDirection;};
  //short  G4int    GetangleIncident() {return angleIncident;};
      
  private:
  
  G4int      CathodId;
  G4double   time;
  G4int      originalInfo;
  G4int      IMany;
#ifdef G4MYLASER_PARAMETERIZATION
  G4int      IManyScatters;
  G4double   ScatteringSteps[100];
  G4double   ScatteringAngles[100];
#endif
  //short  G4int      angleIncident;
  //short  G4int      angleDirection;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<KM3Hit> KM3HitsCollection;

extern G4Allocator<KM3Hit> KM3HitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* KM3Hit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) KM3HitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void KM3Hit::operator delete(void *aHit)
{
  KM3HitAllocator.FreeSingle((KM3Hit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


