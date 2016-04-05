#ifndef EvtIO_h
#define EvtIO_h

#include "io_gcc.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <vector>

using namespace std;

class EvtIO {
 public:
  EvtIO(char *infile, char *outfile);
  ~EvtIO();

 public:
  void ReadRunHeader();
  void WriteRunHeader();
  void ReadEvent();
  void WriteEvent();
  void AddHit(int id, int PMTid, double pe, double t, int trackid, int npepure,
              double ttpure, int creatorProcess);
  void AddNumberOfHits(int hitnumber);
  void AddMuonPositionInfo(int tracknumber, int positionnumber, double posx,
                           double posy, double posz, double momx, double momy,
                           double momz, double mom, double time);
  void AddMuonPositionInfo(int tracknumber, int positionnumber, double posx,
                           double posy, double posz, double time);
  void AddMuonDecaySecondaries(int trackID, int parentID, double posx,
                               double posy, double posz, double dx, double dy,
                               double dz, double energy, double time,
                               int idPDG);
#ifdef G4MYMUON_KEEPENERGY
  void AddMuonEnergyInfo(const vector<double> &info);
#endif

 private:
  event *evt;
  std::ifstream infile;
  std::ofstream outfile;
  bool RunHeaderIsRead;
  bool RunHeaderIsWrite;
  int ParticlesHEPNumber[210000];
  int ParticlesIdNumber[210000];
  bool isneutrinoevent;
  bool hasbundleinfo;
  void GetArgs(string &chd, int &argnumber, double *args);
  int NumberOfParticles;
};
#endif
