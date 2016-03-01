//

#ifndef KM3SD_h
#define KM3SD_h 1

#include "G4VSensitiveDetector.hh"
#include "KM3Hit.hh"
#include <stdio.h>
#include "KM3Detector.hh"
#include "Randomize.hh"
#include "G4MaterialPropertiesTable.hh"

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class KM3SD : public G4VSensitiveDetector {
public:
  FILE *outfile;
  KM3SD(std::string);
  ~KM3SD();

  void Initialize(G4HCofThisEvent *);
  bool ProcessHits(G4Step *, G4TouchableHistory *);
  void EndOfEvent(G4HCofThisEvent *);
  KM3Detector *myStDetector;
  // short  void InsertExternalHit(int ic,double time,int
  // originalInfo,int angleDirection,int angleIncident);
  void InsertExternalHit(int ic, const G4ThreeVector &OMPosition,
                         double time, int originalInfo,
                         const G4ThreeVector &photonDirection);
#if defined(G4MYEM_PARAMETERIZATION) || defined(G4MYHA_PARAMETERIZATION) // newha
  std::vector<long double> *myPhotonsNumber;
  std::vector<long double> *myPhotonsTime;
  std::vector<long double> *myPhotonsTh2Th3Num;
  std::vector<long double> *myPhotonsTh2;
  std::vector<long double> *myPhotonsTh3;
  int NumberOfThetas;
  int NumberOfPhis[47];
  double theta_High[47];
  double phi_High[47][20];
  double theta_Low[47];
  double phi_Low[47][20];
#endif
private:
  KM3HitsCollection *MyCollection;
  int ProcessMyCollection(KM3HitsCollection *aCollection);
  void DrawCathodHit(int NumberOfPhotons, G4ThreeVector pos);
  double TResidual(double, const G4ThreeVector &, const G4ThreeVector &,
                     const G4ThreeVector &);
  void clear();
  void PrintAll();
  void QuickSort(int shorttype, std::vector<KM3Hit *> *theCollectionVector,
                 int top, int bottom);
  int partition_CathodId(std::vector<KM3Hit *> *theCollectionVector,
                           int top, int bottom);
  int partition_Time(std::vector<KM3Hit *> *theCollectionVector, int top,
                       int bottom);
  bool AcceptAngle(double cosangle, double CathodRadius,
                     double CathodHeight, bool);
  void MergeHits(int nfirst, int nlast, double MergeWindow);
  double thespeedmaxQE;
#if defined(G4MYEM_PARAMETERIZATION) || defined(G4MYHA_PARAMETERIZATION) // newha
  G4MaterialPropertyVector *VirtualAbsVector;
  G4MaterialPropertyVector *TrueAbsVector;
#endif
#ifdef G4MYFIT_PARAMETERIZATION
  int ArrayParam[20000];
  int ArrayParamAll[20000];
#endif
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
