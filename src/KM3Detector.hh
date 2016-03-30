#ifndef KM3Detector_h
#define KM3Detector_h 1
#include "KM3Definitions.h"
#include <stdio.h>
#include <vector>
#include "KM3Cathods.hh"
#include "KM3PrimaryGeneratorAction.hh"

#include "HOURSevtWRITE.hh"

class G4VPhysicalVolume;
#include "G4VUserDetectorConstruction.hh"
// newgeant #include "Saxana/SAXProcessor.h"
// newgeant #include "Saxana/ProcessingConfigurator.h"

class KM3Detector : public G4VUserDetectorConstruction {
 public:
  KM3Detector();
  ~KM3Detector();

 public:
  FILE *outfile;
  FILE *outfilePar;
  G4bool useANTARESformat;
  HOURSevtWRITE *TheEVTtoWrite;

  G4VPhysicalVolume *Construct();
  G4double Quantum_Efficiency;
  G4double bottomPosition;
  G4ThreeVector detectorCenter;
  G4double lowestStorey;
  G4double highestStorey;
  G4double outerStorey;
  G4double detectorRadius;

  // this is the maximum vertical distance of the storeys
  // from the center plus a number of absorpion lengths
  G4double detectorMaxz;

  // this is the maximum orizontal distance of the
  // storeys plus a number of absorption lengths
  G4double detectorMaxRho;

  std::vector<StoreysPositions *> *allStoreys;
  std::vector<OMPositions *> *allOMs;
  std::vector<TowersPositions *> *allTowers;  // new towers
  KM3Cathods *allCathods;
  G4double MaxAbsDist;
  G4bool vrmlhits;
  G4bool DrawDetector;
  char *Geometry_File;
  char *Parameter_File;
  G4double TotCathodArea;
  KM3PrimaryGeneratorAction *MyGenerator;

 private:
  void FindDetectorRadius(void);
  void ConstructMaterials(void);
  G4int TotalPMTEntities(const G4VPhysicalVolume *) const;
  void SetUpVariables(void);
  void ReadPMTPositions(const std::string &detx_file_name);
  // newgeant  void sxpInitialize(void);

 private:
  // newgeant  SAXProcessor sxp;
  // newgeant ProcessingConfigurator config;
  G4VPhysicalVolume *fWorld;

 private:
  G4double detectorDepth;
  G4int NUMENTRIES;
  G4int NUMENTRIES_ANGLEACC;
  G4double PPCKOV[100];
  G4double RINDEX_WATER[100];
  G4double Water_Transparency;
  G4double ABSORPTION_WATER[100];
  G4double ABSORPTION_GLASS[100];
  G4double RINDEX_GLASS[100];
  G4double ABSORPTION_GELL[100];
  G4double RINDEX_GELL[100];
  G4double ABSORPTION_AIR[100];
  G4double RINDEX_AIR[100];
  G4double ABSORPTION_CATH[100];
  G4double RINDEX_CATH[100];
  G4double Q_EFF[100];
  G4double SCATTER_WATER[100];
  G4double MieModel;
  G4double COSANGLES[100];
  G4double ACCEPTANCE[100];

  int global_det_id_;
  int n_doms_;

}

};
#endif  // KM3Detector_h
