#include "HOURSevtWRITE.h"

HOURSevtWRITE::HOURSevtWRITE(char *infilechar, char *outfilechar) {
  infile.open(infilechar, ifstream::in);
  evt = new event();

  // the following is to find if it is neutrino events
  int ierr = evt->read(infile);
  int nevents = 0;
  isneutrinoevent = true;
  hasbundleinfo = true;
  while (evt->read(infile) == 0) {
    nevents++;
    if (nevents < 10 && evt->ndat("neutrino") == 0) {
      isneutrinoevent = false;
      break;
    }
    if (nevents < 10 && evt->ndat("track_bundle") == 0) {
      hasbundleinfo = false;
      break;
    }
    if (nevents == 10) break;
  }
  // position to the beggining of the file
  infile.clear();
  infile.seekg(0, ios::beg);
  ////////////////////////////////////////////////////

  outfile.open(outfilechar, ofstream::out);
  RunHeaderIsRead = false;
  RunHeaderIsWrite = false;
}

HOURSevtWRITE::~HOURSevtWRITE() {
  delete evt;
  infile.close();
  outfile.close();
}

void HOURSevtWRITE::ReadRunHeader() {
  if (!RunHeaderIsRead) evt->read(infile);
  RunHeaderIsRead = true;
}

void HOURSevtWRITE::WriteRunHeader() {
  if (!RunHeaderIsWrite) evt->write(outfile);
  RunHeaderIsWrite = true;
}

void HOURSevtWRITE::ReadEvent() {
  evt->read(infile);
  double args[100];
  int argnumber;
  bool UseEarthLepton = false;
  // use earthlepton or not
  if (isneutrinoevent && !hasbundleinfo) {
    evt->ndat("neutrino");
    string NeutrinoInfo = evt->next("neutrino");
    GetArgs(NeutrinoInfo, argnumber, args);
    double xneu, yneu, zneu;
    xneu = args[1];
    yneu = args[2];
    zneu = args[3];
    int NumberOfPart = evt->ndat("track_in");
    if (NumberOfPart > 1) {
      double xx0, yy0, zz0;
      for (int ipart = 0; ipart < NumberOfPart; ipart++) {
        string ParticleInfo = evt->next("track_in");
        GetArgs(ParticleInfo, argnumber, args);
        xx0 = args[1];
        yy0 = args[2];
        zz0 = args[3];
        if (xx0 != xneu || yy0 != yneu || zz0 != zneu) {
          UseEarthLepton = true;
          break;
        }
      }
    }
  }
  if (UseEarthLepton)
    NumberOfParticles = evt->ndat("track_earthlepton");
  else
    NumberOfParticles = evt->ndat("track_in");
  int icount = 0;
  for (int ip = 0; ip < NumberOfParticles; ip++) {
    string ParticleInfo;
    if (UseEarthLepton)
      ParticleInfo = evt->next("track_earthlepton");
    else
      ParticleInfo = evt->next("track_in");
    GetArgs(ParticleInfo, argnumber, args);
    if ((int)args[9] >
        0) {  // do not count particles not within standard pdg coding
      ParticlesIdNumber[icount] = (int)args[0];
      ParticlesHEPNumber[icount] = (int)args[9];
      icount++;
    }
  }
}

void HOURSevtWRITE::GetArgs(string &chd, int &argnumber, double *args) {
  string subchd = chd;
  size_t length = subchd.length();
  size_t start, stop;
  argnumber = 0;
  while (length > 0) {
    start = 0;
    stop = subchd.find_first_of(" ");
    if (stop != string::npos) {
      args[argnumber] = atof((subchd.substr(start, stop - start)).data());
      start = subchd.find_first_not_of(" ", stop);
      if (start != string::npos) {
        subchd = subchd.substr(start, length);
        length = subchd.length();
      } else {
        length = 0;
      }
    } else {
      args[argnumber] = atof(subchd.data());
      length = 0;
    }
    argnumber++;
  }
}

void HOURSevtWRITE::WriteEvent() { evt->write(outfile); }

void HOURSevtWRITE::AddHit(int id, int PMTid, double pe, double t, int trackid,
                           int npepure, double ttpure, int creatorProcess) {
  string dt("hit");
  char buffer[256];
  int Gid;
  if (trackid <= NumberOfParticles) {
    Gid = ParticlesHEPNumber[trackid - 1];
    // convert from geant track id to input track id
    trackid = ParticlesIdNumber[trackid - 1];
  } else {
    Gid = 6;
    trackid = 999999;
  }
  PMTid++;  // in the evt file the numbering of pmts starts from 1
  sprintf(buffer, "%8d %6d %6.2f %10.2f %4d %4d %3d %10.2f %4d", id, PMTid, pe,
          t, Gid, trackid, npepure, ttpure, creatorProcess);
  string dw(buffer);
  evt->taga(dt, dw);
}

void HOURSevtWRITE::AddNumberOfHits(int hitnumber) {
  string dt("total_hits");
  char buffer[256];
  sprintf(buffer, "%8d", hitnumber);
  string dw(buffer);
  evt->taga(dt, dw);
}

void HOURSevtWRITE::AddMuonPositionInfo(int tracknumber, int positionnumber,
                                        double posx, double posy, double posz,
                                        double momx, double momy, double momz,
                                        double mom, double time) {
  string dt("muonaddi_info");
  char buffer[256];
  tracknumber =
      ParticlesIdNumber[tracknumber -
                        1];  // convert from geant track id to input track id
  sprintf(buffer,
          "%4d %2d %8.2f %8.2f %8.2f %10.6f %10.6f %10.6f %12.6e %10.2f",
          tracknumber, positionnumber, posx, posy, posz, momx, momy, momz, mom,
          time);
  string dw(buffer);
  evt->taga(dt, dw);
}

void HOURSevtWRITE::AddMuonPositionInfo(int tracknumber, int positionnumber,
                                        double posx, double posy, double posz,
                                        double time) {
  string dt("muonaddi_info");
  char buffer[256];
  tracknumber =
      ParticlesIdNumber[tracknumber -
                        1];  // convert from geant track id to input track id
  sprintf(buffer, "%4d %2d %8.2f %8.2f %8.2f %10.2f", tracknumber,
          positionnumber, posx, posy, posz, time);
  string dw(buffer);
  evt->taga(dt, dw);
}

void HOURSevtWRITE::AddMuonDecaySecondaries(int trackID, int parentID,
                                            double posx, double posy,
                                            double posz, double dx, double dy,
                                            double dz, double energy,
                                            double time, int idPDG) {
  if (parentID > NumberOfParticles) return;
  int Gid = ParticlesHEPNumber[parentID - 1];
  if ((Gid != 5) && (Gid != 6)) return;
  parentID =
      ParticlesIdNumber[parentID -
                        1];  // convert from geant track id to input track id
  string dt("muon_decay");
  char buffer[256];
  sprintf(
      buffer,
      "%6d %6d %10.3f %10.3f %10.3f %12.8f %12.8f %12.8f %12.6f %10.2f %10d",
      trackID, parentID, posx, posy, posz, dx, dy, dz, energy, time, idPDG);
  string dw(buffer);
  evt->taga(dt, dw);
}

#ifdef G4MYMUON_KEEPENERGY
void HOURSevtWRITE::AddMuonEnergyInfo(const vector<double> &info) {
  if (info.size() == 0) return;
  string dt("muonenergy_info");
  char buffer[256];
  int NumberOfTags = int((info.size() - 1) / 10.0) + 1;
  for (int itag = 0; itag < NumberOfTags; itag++) {
    int istart = 10 * itag;
    int istop = istart + 10;
    if (istop > info.size()) istop = info.size();
    int nentries = istop - istart;
    int i = istart;
    if (nentries == 10)
      sprintf(buffer,
              "%4d %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e",
              itag, info[i], info[i + 1], info[i + 2], info[i + 3], info[i + 4],
              info[i + 5], info[i + 6], info[i + 7], info[i + 8], info[i + 9]);
    else if (nentries == 9)
      sprintf(buffer,
              "%4d %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e", itag,
              info[i], info[i + 1], info[i + 2], info[i + 3], info[i + 4],
              info[i + 5], info[i + 6], info[i + 7], info[i + 8]);
    else if (nentries == 8)
      sprintf(buffer, "%4d %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e",
              itag, info[i], info[i + 1], info[i + 2], info[i + 3], info[i + 4],
              info[i + 5], info[i + 6], info[i + 7]);
    else if (nentries == 7)
      sprintf(buffer, "%4d %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e", itag,
              info[i], info[i + 1], info[i + 2], info[i + 3], info[i + 4],
              info[i + 5], info[i + 6]);
    else if (nentries == 6)
      sprintf(buffer, "%4d %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e", itag, info[i],
              info[i + 1], info[i + 2], info[i + 3], info[i + 4], info[i + 5]);
    else if (nentries == 5)
      sprintf(buffer, "%4d %8.2e %8.2e %8.2e %8.2e %8.2e", itag, info[i],
              info[i + 1], info[i + 2], info[i + 3], info[i + 4]);
    else if (nentries == 4)
      sprintf(buffer, "%4d %8.2e %8.2e %8.2e %8.2e", itag, info[i], info[i + 1],
              info[i + 2], info[i + 3]);
    else if (nentries == 3)
      sprintf(buffer, "%4d %8.2e %8.2e %8.2e", itag, info[i], info[i + 1],
              info[i + 2]);
    else if (nentries == 2)
      sprintf(buffer, "%4d %8.2e %8.2e", itag, info[i], info[i + 1]);
    else if (nentries == 1)
      sprintf(buffer, "%4d %8.2e", itag, info[i]);
    string dw(buffer);
    evt->taga(dt, dw);
  }
}
#endif
