#ifndef HOURSevtREAD_h
#define HOURSevtREAD_h

// IS DEPRECATED
// IO IS HAPPENING IN
// EvtIO

#include "io_gcc.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
// following is for hepevt interface
#include "G4Event.h"
//////////////////////////////////
using namespace std;

class HOURSevtREAD {
 public:
  HOURSevtREAD(char *infile);
  ~HOURSevtREAD();

 public:
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
  event *evt;
  int nevents;
  ifstream infile;
  int ICONPDG[174];
  double PDGMASS[174];
  void Initialize(void);
  int ConvertHEPToPDG(int hepcode);
  void GetArgs(string &chd, int &argnumber, double *args);
  double GetParticleMass(int hepcode);
  bool UseEarthLepton;
  bool isneutrinoevent;
  bool hasbundleinfo;
  bool ReadNeutrinoVertexParticles;
};
#endif
