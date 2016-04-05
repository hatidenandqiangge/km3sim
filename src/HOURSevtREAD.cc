#include "HOURSevtREAD.h"

HOURSevtREAD::HOURSevtREAD(char *infilechar) {
  infile.open(infilechar, std::ifstream::in);
  evt = new event();
  // read header
  int ierr = evt->read(infile);
  // count events
  nevents = 0;
  isneutrinoevent = true;
  hasbundleinfo = true;
  while (evt->read(infile) == 0) {
    nevents++;
    if (nevents < 10 && evt->ndat("neutrino") == 0) isneutrinoevent = false;
    if (nevents < 10 && evt->ndat("track_bundle") == 0) hasbundleinfo = false;
  }
  // position to the beggining of the file
  infile.clear();
  infile.seekg(0, std::ios::beg);

  // read header again
  ierr = evt->read(infile);

  Initialize();
}

HOURSevtREAD::~HOURSevtREAD() {
  delete evt;
  infile.close();
}

int HOURSevtREAD::GetNumberOfEvents() { return nevents; }

void HOURSevtREAD::ReadEvent(void) {
  evt->read(infile);
  UseEarthLepton = false;
  if (isneutrinoevent && !hasbundleinfo) {
    int idneu, idtarget;
    double xneu, yneu, zneu, pxneu, pyneu, pzneu, t0;
    GetNeutrinoInfo(idneu, idtarget, xneu, yneu, zneu, pxneu, pyneu, pzneu, t0);
    int NumberOfParticles = GetNumberOfParticles();
    if (NumberOfParticles > 1) {
      int idbeam;
      double xx0, yy0, zz0, pxx0, pyy0, pzz0, t0;
      for (int ipart = 0; ipart < NumberOfParticles; ipart++) {
        GetParticleInfo(idbeam, xx0, yy0, zz0, pxx0, pyy0, pzz0, t0);
        if (xx0 != xneu || yy0 != yneu || zz0 != zneu) {
          UseEarthLepton = true;
          break;
        }
      }
    }
  }
}

int HOURSevtREAD::GetNumberOfParticles(void) {
  if (UseEarthLepton)
    return evt->ndat("track_earthlepton");
  else
    return evt->ndat("track_in");
}

void HOURSevtREAD::Initialize(void) {
  // convert hep to pdg
  for (int i = 0; i <= 173; i++) ICONPDG[i] = 0;

  ICONPDG[1] = 22;
  ICONPDG[2] = -11;
  ICONPDG[3] = 11;

  ICONPDG[5] = -13;
  ICONPDG[6] = 13;
  ICONPDG[7] = 111;
  ICONPDG[8] = 211;
  ICONPDG[9] = -211;
  ICONPDG[10] = 130;
  ICONPDG[11] = 321;
  ICONPDG[12] = -321;
  ICONPDG[13] = 2112;
  ICONPDG[14] = 2212;
  ICONPDG[15] = -2212;
  ICONPDG[16] = 310;
  ICONPDG[17] = 221;
  ICONPDG[18] = 3122;
  ICONPDG[19] = 3222;
  ICONPDG[20] = 3212;
  ICONPDG[21] = 3112;
  ICONPDG[22] = 3322;
  ICONPDG[23] = 3312;
  ICONPDG[24] = 3334;
  ICONPDG[25] = -2112;
  ICONPDG[26] = -3122;
  ICONPDG[27] = -3112;
  ICONPDG[28] = -3212;
  ICONPDG[29] = -3222;
  ICONPDG[30] = -3322;
  ICONPDG[31] = -3312;
  ICONPDG[32] = -3334;

  ICONPDG[50] = 223;
  ICONPDG[51] = 113;
  ICONPDG[52] = 213;
  ICONPDG[53] = -213;
  ICONPDG[54] = 2224;
  ICONPDG[55] = 2214;
  ICONPDG[56] = 2114;
  ICONPDG[57] = 1114;
  ICONPDG[58] = -2224;
  ICONPDG[59] = -2214;
  ICONPDG[60] = -2114;
  ICONPDG[61] = -1114;
  ICONPDG[62] = 313;
  ICONPDG[63] = 323;
  ICONPDG[64] = -323;
  ICONPDG[65] = -313;
  ICONPDG[66] = 12;
  ICONPDG[67] = -12;
  ICONPDG[68] = 14;
  ICONPDG[69] = -14;

  ICONPDG[116] = 421;
  ICONPDG[117] = 411;
  ICONPDG[118] = -411;
  ICONPDG[119] = -421;
  ICONPDG[120] = 431;
  ICONPDG[121] = -431;
  ICONPDG[122] = 441;
  ICONPDG[123] = 423;
  ICONPDG[124] = 413;
  ICONPDG[125] = -413;
  ICONPDG[126] = -423;
  ICONPDG[127] = 433;
  ICONPDG[128] = -433;

  ICONPDG[130] = 443;
  ICONPDG[131] = -15;
  ICONPDG[132] = 15;
  ICONPDG[133] = 16;
  ICONPDG[134] = -16;

  ICONPDG[137] = 4122;
  ICONPDG[138] = 4232;
  ICONPDG[139] = 4132;
  ICONPDG[140] = 4222;
  ICONPDG[141] = 4212;
  ICONPDG[142] = 4112;
  ICONPDG[143] = 4322;
  ICONPDG[144] = 4312;
  ICONPDG[145] = 4332;

  ICONPDG[149] = -4122;
  ICONPDG[150] = -4232;
  ICONPDG[151] = -4132;
  ICONPDG[152] = -4222;
  ICONPDG[153] = -4212;
  ICONPDG[154] = -4112;
  ICONPDG[155] = -4322;
  ICONPDG[156] = -4312;
  ICONPDG[157] = -4332;

  ICONPDG[161] = 4224;
  ICONPDG[162] = 4214;
  ICONPDG[163] = 4114;

  ICONPDG[171] = -4224;
  ICONPDG[172] = -4214;
  ICONPDG[173] = -4114;

  // hep masses
  for (int i = 0; i <= 173; i++) PDGMASS[i] = 0.0;

  PDGMASS[1] = 0.0;
  PDGMASS[2] = 0.510998902E-3;
  PDGMASS[3] = 0.510998902E-3;
  PDGMASS[5] = 0.105658357;
  PDGMASS[6] = 0.105658357;
  PDGMASS[7] = 0.1349766;
  PDGMASS[8] = 0.13957018;
  PDGMASS[9] = 0.13957018;
  PDGMASS[10] = 0.497672;
  PDGMASS[11] = 0.493677;
  PDGMASS[12] = 0.493677;
  PDGMASS[13] = 0.93956533;
  PDGMASS[14] = 0.93827200;
  PDGMASS[15] = 0.93827200;
  PDGMASS[16] = 0.497672;
  PDGMASS[17] = 0.54730;
  PDGMASS[18] = 1.115683;
  PDGMASS[19] = 1.18937;
  PDGMASS[20] = 1.192642;
  PDGMASS[21] = 1.197449;
  PDGMASS[22] = 1.31483;
  PDGMASS[23] = 1.32131;
  PDGMASS[24] = 1.67245;
  PDGMASS[25] = 0.93956533;
  PDGMASS[26] = 1.115683;
  PDGMASS[27] = 1.18937;
  PDGMASS[28] = 1.192642;
  PDGMASS[29] = 1.197449;
  PDGMASS[30] = 1.31483;
  PDGMASS[31] = 1.32131;
  PDGMASS[32] = 1.67245;

  PDGMASS[50] = 0.78257;
  PDGMASS[51] = 0.7690;
  PDGMASS[52] = 0.7665;
  PDGMASS[53] = 0.7665;
  PDGMASS[54] = 1.2305;
  PDGMASS[55] = 1.2318;
  PDGMASS[56] = 1.2331;
  PDGMASS[57] = 1.2344;
  PDGMASS[58] = 1.2309;
  PDGMASS[59] = 1.2323;
  PDGMASS[60] = 1.2336;
  PDGMASS[61] = 1.2349;
  PDGMASS[62] = 0.89610;
  PDGMASS[63] = 0.89166;
  PDGMASS[64] = 0.89166;
  PDGMASS[65] = 0.89610;
  PDGMASS[66] = 0.0;
  PDGMASS[67] = 0.0;
  PDGMASS[68] = 0.0;
  PDGMASS[69] = 0.0;

  PDGMASS[116] = 1.8645;
  PDGMASS[117] = 1.8693;
  PDGMASS[118] = 1.8693;
  PDGMASS[119] = 1.8645;
  PDGMASS[120] = 1.9686;
  PDGMASS[121] = 1.9685;
  PDGMASS[122] = 2.9797;
  PDGMASS[123] = 2.0067;
  PDGMASS[124] = 2.0100;
  PDGMASS[125] = 2.0100;
  PDGMASS[126] = 2.0067;
  PDGMASS[127] = 2.1124;
  PDGMASS[128] = 2.1124;

  PDGMASS[130] = 3.09687;
  PDGMASS[131] = 1.77699;
  PDGMASS[132] = 1.77699;
  PDGMASS[133] = 0.0;
  PDGMASS[134] = 0.0;

  PDGMASS[137] = 2.2849;
  PDGMASS[138] = 2.4663;
  PDGMASS[139] = 2.4718;
  PDGMASS[140] = 2.4526;
  PDGMASS[141] = 2.4513;
  PDGMASS[142] = 2.4522;
  PDGMASS[143] = 2.5741;
  PDGMASS[144] = 2.5788;
  PDGMASS[145] = 2.6975;

  PDGMASS[149] = 2.2849;
  PDGMASS[150] = 2.4663;
  PDGMASS[151] = 2.4718;
  PDGMASS[152] = 2.4526;
  PDGMASS[153] = 2.4513;
  PDGMASS[154] = 2.4522;
  PDGMASS[155] = 2.5741;
  PDGMASS[156] = 2.5788;
  PDGMASS[157] = 2.6975;

  PDGMASS[161] = 2.5194;
  PDGMASS[162] = 2.5159;
  PDGMASS[163] = 2.5175;

  PDGMASS[171] = 2.5194;
  PDGMASS[172] = 2.5159;
  PDGMASS[173] = 2.5175;
}

int HOURSevtREAD::ConvertHEPToPDG(int hepcode) { return ICONPDG[hepcode]; }

double HOURSevtREAD::GetParticleMass(int hepcode) { return PDGMASS[hepcode]; }

void HOURSevtREAD::GetParticleInfo(int &idbeam, double &xx0, double &yy0,
                                   double &zz0, double &pxx0, double &pyy0,
                                   double &pzz0, double &t0) {
  std::string ParticleInfo;
  if (UseEarthLepton)
    ParticleInfo = evt->next("track_earthlepton");
  else
    ParticleInfo = evt->next("track_in");
  double args[100];
  int argnumber;
  GetArgs(ParticleInfo, argnumber, args);
  if ((int)args[9] <= 0) {  // in order to get rid off particles that are not
    // standard (pythia or genie internal code particles
    // e.g. 93)
    idbeam = 0;
    return;
  }
  if (isneutrinoevent && hasbundleinfo) {  // in order to select particles from
    // neutrino interaction or muon bundle
    if (ReadNeutrinoVertexParticles &&
        (int)args[11] == 1) {  // selsct neutrino interaction particles
      idbeam = 0;
      return;
    }
    if (!ReadNeutrinoVertexParticles &&
        (int)args[11] == 0) {  // selsct bundle muons
      idbeam = 0;
      return;
    }
  }
  if (argnumber > 10)
    idbeam = (int)args[10];
  else
    idbeam = ConvertHEPToPDG((int)args[9]);
  xx0 = args[1] * m / cm;  // convert to cm
  yy0 = args[2] * m / cm;
  zz0 = args[3] * m / cm;
  double totenergy = args[7];
  double pmass = GetParticleMass((int)args[9]);
  double pmom = totenergy * totenergy - pmass * pmass;
  if (pmom > 0.0)
    pmom = sqrt(totenergy * totenergy - pmass * pmass);
  else
    pmom = 0.0;
  pxx0 = args[4] * pmom;
  pyy0 = args[5] * pmom;
  pzz0 = args[6] * pmom;
  t0 = args[8];
}

void HOURSevtREAD::GetNeutrinoInfo(int &idneu, int &idtarget, double &xneu,
                                   double &yneu, double &zneu, double &pxneu,
                                   double &pyneu, double &pzneu, double &t0) {
  idneu = 0;
  idtarget = 0;
  xneu = 0.0;
  yneu = 0.0;
  zneu = 0.0;
  pxneu = 0.0;
  pyneu = 0.0;
  pzneu = 0.0;
  if (!isneutrinoevent) return;
  evt->ndat("neutrino");
  std::string NeutrinoInfo = evt->next("neutrino");
  double args[100];
  int argnumber;
  GetArgs(NeutrinoInfo, argnumber, args);
  idneu = (int)args[12];
  if (argnumber == 15)
    idtarget = (int)args[14];
  else
    idtarget = 1445;
  xneu = args[1] * m / cm;  // convert to cm
  yneu = args[2] * m / cm;
  zneu = args[3] * m / cm;
  double pmom = args[7];
  pxneu = args[4] * pmom;
  pyneu = args[5] * pmom;
  pzneu = args[6] * pmom;
  t0 = args[8];
  // in case that ther is one track_in not with the same vertex with neutrino
  if (evt->ndat("track_in") == 1 && !hasbundleinfo) {
    int idbeam;
    double xx0, yy0, zz0;
    double pxx0, pyy0, pzz0;
    double t0;
    GetParticleInfo(idbeam, xx0, yy0, zz0, pxx0, pyy0, pzz0, t0);
    if (xx0 != xneu || yy0 != yneu || zz0 != zneu) {
      xneu = xx0;
      yneu = yy0;
      zneu = zz0;
    }
  }
}

void HOURSevtREAD::GetArgs(std::string &chd, int &argnumber, double *args) {
  std::string subchd = chd;
  size_t length = subchd.length();
  size_t start, stop;
  argnumber = 0;
  while (length > 0) {
    start = 0;
    stop = subchd.find_first_of(" ");
    if (stop != std::string::npos) {
      args[argnumber] = atof((subchd.substr(start, stop - start)).data());
      start = subchd.find_first_not_of(" ", stop);
      if (start != std::string::npos) {
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

bool HOURSevtREAD::IsNeutrinoEvent(void) { return isneutrinoevent; }

// following is for hepevt interface
#include "G4PrimaryVertex.h"
#include "G4PrimaryParticle.h"
//#include "G4ThreeVector.h"

//////////////////////////////////

void HOURSevtREAD::GeneratePrimaryVertex(G4Event *anEvent) {
  if (isneutrinoevent && hasbundleinfo) {
    // first read the information of the neutrino
    // vertex////////////////////////////////
    int idneu, idtarget;
    double xneu, yneu, zneu, pxneu, pyneu, pzneu, t0;
    GetNeutrinoInfo(idneu, idtarget, xneu, yneu, zneu, pxneu, pyneu, pzneu, t0);
    // create the neutrino interaction vertex including only the time
    // the position is set in generator action
    G4ThreeVector zero;
    G4PrimaryVertex *vertex = new G4PrimaryVertex(zero, t0);
    // add the particles from neutrino interaction to the vertex
    int NHEP;  // number of entries
    NHEP = GetNumberOfParticles();
    ReadNeutrinoVertexParticles = true;
    int idbeam;
    double xx0, yy0, zz0;
    double pxx0, pyy0, pzz0;
    for (int IHEP = 0; IHEP < NHEP; IHEP++) {
      GetParticleInfo(idbeam, xx0, yy0, zz0, pxx0, pyy0, pzz0, t0);
      if (idbeam != 0) {  // do not load particles not within PDG coding (e.g.
        // 93) and from muon bundle
        if (abs(idbeam) != 411 && abs(idbeam) != 421 && abs(idbeam) != 431 &&
            abs(idbeam) != 4122 && abs(idbeam) != 4212 &&
            abs(idbeam) != 4222) {  // these particles are not defined or have
          // not decay modes in GEANT4
          G4PrimaryParticle *particle = new G4PrimaryParticle(idbeam);
          particle->SetMomentum(pxx0 * GeV, pyy0 * GeV, pzz0 * GeV);
          vertex->SetPrimary(particle);
        }
      }
    }
    // Put the vertex to G4Event object
    anEvent->AddPrimaryVertex(vertex);
    //////////////////////////////////////////////////////////////////////////////////
    // next load the information from the bundle
    // muons/////////////////////////////////
    NHEP = GetNumberOfParticles();
    ReadNeutrinoVertexParticles = false;
    for (int IHEP = 0; IHEP < NHEP; IHEP++) {
      GetParticleInfo(idbeam, xx0, yy0, zz0, pxx0, pyy0, pzz0, t0);
      if (idbeam != 0) {  // load particles from muon bundle only
        G4PrimaryParticle *particle = new G4PrimaryParticle(idbeam);
        particle->SetMomentum(pxx0 * GeV, pyy0 * GeV, pzz0 * GeV);
        vertex = new G4PrimaryVertex(
            G4ThreeVector(xx0 * cm, yy0 * cm, zz0 * cm), t0 * ns);
        vertex->SetPrimary(particle);
        anEvent->AddPrimaryVertex(vertex);
      }
    }
  } else {
    // create G4PrimaryVertex object
    G4ThreeVector zero;
    double particle_time = 0.0;
    G4PrimaryVertex *vertex = new G4PrimaryVertex(zero, particle_time);

    int NHEP;  // number of entries
    NHEP = GetNumberOfParticles();
    for (int IHEP = 0; IHEP < NHEP; IHEP++) {
      int idbeam;
      double xx0, yy0, zz0;
      double pxx0, pyy0, pzz0;
      double t0;
      GetParticleInfo(idbeam, xx0, yy0, zz0, pxx0, pyy0, pzz0, t0);
      if (idbeam != 0) {  // do not load particles not within PDG coding (e.g.
        // 93)
        if (abs(idbeam) != 411 && abs(idbeam) != 421 && abs(idbeam) != 431 &&
            abs(idbeam) != 4122 && abs(idbeam) != 4212 &&
            abs(idbeam) != 4222) {  // these particles are not defined or have
          // not decay modes in GEANT4
          G4PrimaryParticle *particle = new G4PrimaryParticle(idbeam);
          particle->SetMomentum(pxx0 * GeV, pyy0 * GeV, pzz0 * GeV);
          vertex->SetPrimary(particle);
        }
      }
    }
    // Put the vertex to G4Event object
    anEvent->AddPrimaryVertex(vertex);
  }
}
