#ifndef HOURSevtREAD_h
#define HOURSevtREAD_h

#include <math.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// following is for hepevt interface
#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
//#include "G4ThreeVector.h"

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include "seaweed.h"

class HOURSevtRead {
 public:
  //HOURSevtRead(char *infile);
  HOURSevtRead(std::string infile);
  ~HOURSevtRead();

  int GetNumberOfEvents();
  void ReadEvent(void);
  void GetNeutrinoInfo(int &idneu, int &idtarget, double &xneu, double &yneu,
                       double &zneu, double &pxneu, double &pyneu,
                       double &pzneu, double &t0);
  int GetNumberOfParticles();
  void GetParticleInfo(int &idbeam, double &xx0, double &yy0, double &zz0,
                       double &pxx0, double &pyy0, double &pzz0, double &t0);
  bool IsNeutrinoEvent(void);

  void GeneratePrimaryVertex(G4Event *anEvent);

 private:
  seaweed::event *evt;
  int nevents;
  std::ifstream infile;
  int ICONPDG[174];
  double PDGMASS[174];
  void InitPDGTables(void);
  int ConvertHEPToPDG(int hepcode);
  void GetArgs(std::string &chd, int &argnumber, double *args);
  double GetParticleMass(int hepcode);
  bool UseEarthLepton;
  bool isneutrinoevent;
  bool hasbundleinfo;
  bool ReadNeutrinoVertexParticles;
};
#endif
