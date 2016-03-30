#include "KM3Sim.h"
#include "KM3Physics.hh"
#include "KM3PrimaryGeneratorAction.hh"
#include "KM3StackingAction.hh"
#include "KM3TrackingAction.hh"
#include "KM3SteppingAction.hh"
#include "KM3EventAction.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4UnitsTable.hh"
#ifdef G4VIS_USE
// ftest #include "KM3Visualisation.hh"
#include "G4VisExecutive.hh"
#include "G4TrajectoryDrawByParticleID.hh"
#include "G4TrajectoryParticleFilter.hh"
#endif
#include <stdio.h>

#ifdef G4MYHAMUONS_PARAMETERIZATION
#include <iostream>
#include <fstream>
#include <iomanip>
#endif

int main(int argc, char *argv[]) {
  G4long myseed = atol(argv[2]);
  G4cout << myseed << G4endl;
  CLHEP::HepRandom::setTheSeed(myseed);

  char *Geometry_File = argv[4];
  char *Parameter_File = argv[5];
  char *EMParametrization_FILE = argv[6];
  char *HAParametrization_FILE = argv[7];
  G4bool useHEPEvt;
  G4bool useANTARESformat;
  char *fileParticles;
  char *filePythiaParticles;
  char *fileParamHAmuons = NULL;
  FILE *outfilePar;
  G4double ParamEnergy;
  G4int ParamNumber;
  G4int ParamParticle;
  if ((argv[8][0] == 'P') && (argv[8][1] == 'a') && (argv[8][2] == 'r') &&
      (argv[8][3] == 'a') && (argv[8][4] == 'm') && (argv[8][5] == 'H') &&
      (argv[8][6] == 'A') && (argv[8][7] == '\0')) {
    useHEPEvt = false;
    useANTARESformat = false;
    if ((outfilePar = fopen(argv[9], "w")) == NULL) {
      printf("Error open output Param file\n");
      return 1;
    }
  }
  ParamEnergy = GeV * atof(argv[10]);
  ParamNumber = atoi(argv[11]);
  ParamParticle = atoi(argv[12]);
}
else if ((argv[8][0] == 'P') && (argv[8][1] == 'a') && (argv[8][2] == 'r') &&
         (argv[8][3] == 'a') && (argv[8][4] == 'm') && (argv[8][5] == 'H') &&
         (argv[8][6] == 'A') && (argv[8][7] == 'M') && (argv[8][8] == 'u') &&
         (argv[8][9] == 'o') && (argv[8][10] == 'n') && (argv[8][11] == 's') &&
         (argv[8][12] == '\0')) {
  useHEPEvt = true;
  useANTARESformat = false;
  fileParticles = argv[9];
  filePythiaParticles = argv[10];
  fileParamHAmuons = argv[11];  // this is fed to stacking and generator
  if ((outfilePar = fopen(argv[12], "w")) ==
      NULL) {  // this works the same way as in ParamEM
    printf(
        "Error open output Param photon file for Hadronic "
        "parameterization\n");
    return 1;
  }
}
else if ((argv[8][0] == 'P') && (argv[8][1] == 'a') && (argv[8][2] == 'r') &&
         (argv[8][3] == 'a') && (argv[8][4] == 'm') && (argv[8][5] == 'F') &&
         (argv[8][6] == 'i') && (argv[8][7] == 't') && (argv[8][8] == '\0')) {
  useHEPEvt = false;
  useANTARESformat = false;
  if (argv[9] == NULL) {
    G4cout << "You must give an output param file" << G4endl;
    return 1;
  } else {
    if ((outfilePar = fopen(argv[9], "w")) == NULL) {
      printf("Error open output Param file\n");
      return 1;
    }
  }
  ParamEnergy = GeV * atof(argv[10]);
  ParamNumber = atoi(argv[11]);
}
else if ((argv[8][0] == 'P') && (argv[8][1] == 'a') && (argv[8][2] == 'r') &&
         (argv[8][3] == 'a') && (argv[8][4] == 'm') && (argv[8][5] == 'E') &&
         (argv[8][6] == 'M') && (argv[8][7] == '\0')) {
  useHEPEvt = false;
  useANTARESformat = false;
  if ((outfilePar = fopen(argv[9], "w")) == NULL) {
    printf("Error open output Param file\n");
    return 1;
  }
  ParamEnergy = GeV * atof(argv[10]);
  ParamNumber = atoi(argv[11]);
  ParamParticle = atoi(argv[12]);
}
else if ((argv[8][0] == 'A') && (argv[8][1] == 'N') && (argv[8][2] == 'T') &&
         (argv[8][3] == 'A') && (argv[8][4] == 'R') && (argv[8][5] == 'E') &&
         (argv[8][6] == 'S') && (argv[8][7] == '_') && (argv[8][8] == 'E') &&
         (argv[8][9] == 'V') && (argv[8][10] == 'T') && (argv[8][11] == '_') &&
         (argv[8][12] == 'F') && (argv[8][13] == 'O') && (argv[8][14] == 'R') &&
         (argv[8][15] == 'M') && (argv[8][16] == 'A') && (argv[8][17] == 'T') &&
         (argv[8][18] == '\0')) {
  useHEPEvt = true;
  useANTARESformat = true;
  fileParticles = argv[9];
}
else {
  useHEPEvt = false;
  useANTARESformat = false;
  fileParticles = argv[8];
}
#ifdef G4MYHAMUONS_PARAMETERIZATION
if (fileParamHAmuons == NULL) {
  G4cout << "G4MYHA_PARAMETERIZATION preprocessor flag must be combined with "
            "ParamHA 7th command line argument"
         << G4endl;
  return 1;
}
#endif

FILE *savefile;
HOURSevtWRITE *TheEVTtoWrite;
if (!useANTARESformat) {
  if ((savefile = fopen(argv[3], "w")) == NULL) {
    printf("Error open file\n");
    return 1;
  }
} else {
  TheEVTtoWrite = new HOURSevtWRITE(fileParticles, argv[3]);
}

G4RunManager *runManager = new G4RunManager;

KM3Detector *Mydet = new KM3Detector;
runManager->SetUserInitialization(Mydet);
KM3Physics *MyPhys = new KM3Physics;
runManager->SetUserInitialization(MyPhys);
MyPhys->aDetector = Mydet;
Mydet->vrmlhits = visuale;
Mydet->DrawDetector = visual;
Mydet->Geometry_File = Geometry_File;
Mydet->Parameter_File = Parameter_File;
Mydet->EMParametrization_FILE = EMParametrization_FILE;
Mydet->HAParametrization_FILE = HAParametrization_FILE;
Mydet->outfilePar = outfilePar;
#if defined(G4MYEM_PARAMETERIZATION) || \
    defined(G4MYHA_PARAMETERIZATION)  // newha
std::vector<long double> *myPhotonsTime = new std::vector<long double>;
std::vector<long double> *myPhotonsNumber = new std::vector<long double>;
std::vector<long double> *myCumPhotons = new std::vector<long double>;
std::vector<long double> *myCumPhotonsRms = new std::vector<long double>;
std::vector<long double> *myCumNorma = new std::vector<long double>;
std::vector<long double> *myPhotonsTh2Th3Num = new std::vector<long double>;
std::vector<long double> *myPhotonsTh2 = new std::vector<long double>;
std::vector<long double> *myPhotonsTh3 = new std::vector<long double>;
G4bool FineBin = false;
G4int VertexDistanceBins = 40;
G4int VertexSolidAngleBins = 51;
G4int TimeSolidAngleBins = 52;
G4int TimeBins = 111;
G4int OMSolidAngleBins = 834;  // for 834 see KM3SD
if (ParamEnergy == 0.0) {
  FineBin = true;
  VertexSolidAngleBins = 71;
}
G4int VertexDistAngleBins = VertexDistanceBins * VertexSolidAngleBins;
long double zero = 0.0;
for (G4int i = 0; i < VertexDistAngleBins; i++) {
  myCumPhotons->push_back(zero);
  myCumPhotonsRms->push_back(zero);
  myCumNorma->push_back(zero);
  myPhotonsNumber->push_back(zero);
  for (G4int k = 0; k < TimeSolidAngleBins; k++) {
    for (G4int j = 0; j < TimeBins; j++) myPhotonsTime->push_back(zero);
  }
  for (G4int k = 0; k < OMSolidAngleBins; k++) {
    myPhotonsTh2Th3Num->push_back(zero);
    myPhotonsTh2->push_back(zero);
    myPhotonsTh3->push_back(zero);
  }
}
Mydet->myPhotonsTime = myPhotonsTime;
Mydet->myPhotonsNumber = myPhotonsNumber;
Mydet->myPhotonsTh2Th3Num = myPhotonsTh2Th3Num;
Mydet->myPhotonsTh2 = myPhotonsTh2;
Mydet->myPhotonsTh3 = myPhotonsTh3;
#endif

#ifdef G4VIS_USE
// Visualization, if you choose to have it!
// ftest  G4VisManager* visManager = new KM3Visualisation;
// ftest  visManager->Initialize();

G4VisManager *visManager = new G4VisExecutive;
visManager->Initialize();
G4TrajectoryDrawByParticleID *model =
    new G4TrajectoryDrawByParticleID("KM3TrajectoryModel");
model->SetDefault(G4Colour(1.0, 1.0, 1.0, 0.0));  // transparent white
model->Set("mu-", "red");                         // assign red to mu-
model->Set("mu+", "red");                         // also to mu+
//  model->Set("e+", "magenta");
// model->Set("e-", "magenta");
//  model->Set("gamma", "green");
// model->Set("opticalphoton", "green");
visManager->RegisterModel(model);
visManager->SelectTrajectoryModel(model->Name());
G4TrajectoryParticleFilter *modelFilter =
    new G4TrajectoryParticleFilter("KM3TrajectoryFilterModel");
modelFilter->Add("mu-");  // show the mu- trajectories
modelFilter->Add("mu+");  // also the mu+
//  modelFilter->Add("e+");
//  modelFilter->Add("e-");
//  modelFilter->Add("gamma");
// modelFilter->Add("opticalphoton");
visManager->RegisterModel(modelFilter);
#endif

// set mandatory user action class
runManager->SetNumberOfEventsToBeStored(0);
// myGeneratorAction and MyTrackingAction
KM3PrimaryGeneratorAction *myGeneratorAction = new KM3PrimaryGeneratorAction;
myGeneratorAction->outfile = savefile;
myGeneratorAction->useHEPEvt = useHEPEvt;
myGeneratorAction->useANTARESformat = useANTARESformat;
myGeneratorAction->fileParticles = fileParticles;
myGeneratorAction->filePythiaParticles = filePythiaParticles;
myGeneratorAction->ParamEnergy = ParamEnergy;
myGeneratorAction->idbeam = ParamParticle;
Mydet->MyGenerator = myGeneratorAction;

KM3TrackingAction *myTracking = new KM3TrackingAction;
myTracking->TheEVTtoWrite = TheEVTtoWrite;
myTracking->useANTARESformat = useANTARESformat;
myGeneratorAction->myTracking =
    myTracking;  // link between generator and tracking (to provide number of
// initial particles to trackingAction
myGeneratorAction->Initialize();
runManager->SetUserAction(myGeneratorAction);

KM3EventAction *event_action = new KM3EventAction;
runManager->SetUserAction(event_action);
event_action->outfile = savefile;
event_action->TheEVTtoWrite = TheEVTtoWrite;
event_action->useANTARESformat = useANTARESformat;
#if !defined(G4MYEM_PARAMETERIZATION) && \
    !defined(G4MYHA_PARAMETERIZATION)            // newha
myGeneratorAction->event_action = event_action;  // generator knows event to
// set the number of initial
// particles
#endif
#ifdef G4MYFIT_PARAMETERIZATION
Mydet->event_action = event_action;  // KM3SD must know the momentum of the
// particle when crosses the center of the
// detector
#endif
#if defined(G4MYEM_PARAMETERIZATION) || \
    defined(G4MYHA_PARAMETERIZATION)  // newha
event_action->MyGenerator =
    myGeneratorAction;  // event action knows generator for EM param purposes
event_action->MyStDetector = Mydet;
event_action->myPhotonsNumber = myPhotonsNumber;
event_action->myCumPhotons = myCumPhotons;
event_action->myCumPhotonsRms = myCumPhotonsRms;
event_action->myCumNorma = myCumNorma;
#endif

Mydet->outfile = savefile;
Mydet->TheEVTtoWrite = TheEVTtoWrite;
Mydet->useANTARESformat = useANTARESformat;

KM3StackingAction *myStacking = new KM3StackingAction;
KM3SteppingAction *myStepping = new KM3SteppingAction;
myStacking->SetDetector(Mydet);
#ifdef G4MYHAMUONS_PARAMETERIZATION
std::ofstream *outMuonHAFile =
    new std::ofstream(fileParamHAmuons, std::ios::out | std::ios::binary);
myGeneratorAction->outMuonHAFile = outMuonHAFile;
myStacking->outMuonHAFile = outMuonHAFile;
#endif
myStepping->myStDetector = Mydet;
myStepping->event_action = event_action;
runManager->SetUserAction(myStacking);
runManager->SetUserAction(myTracking);
runManager->SetUserAction(myStepping);
// Initialize G4 kernel
runManager->Initialize();
//   for (G4int iom=0 ; iom < Mydet->allOMs->size() ; iom++){
//     G4cout << G4BestUnit((*(Mydet->allOMs))[iom]->position(0),"Length")
//     << G4BestUnit((*(Mydet->allOMs))[iom]->position(1),"Length")
//     << G4BestUnit((*(Mydet->allOMs))[iom]->position(2),"Length")
//     << G4BestUnit((*(Mydet->allOMs))[iom]->radius,"Length") <<"   "<< iom
// <<"\n";
//   }

// get the pointer to the UI manager and set verbosities
G4UImanager *UI = G4UImanager::GetUIpointer();
G4UIsession *session = 0;
#ifdef G4UI_USE_TCSH
session = new G4UIterminal(new G4UItcsh);
#else
session = new G4UIterminal();
#endif

// UI->ApplyCommand("/control/execute myrun.mac");

//  UI->ApplyCommand("/event/verbose 0");
//  UI->ApplyCommand("/control/verbose 0");
//    UI->ApplyCommand("/tracking/verbose 1");
//    UI->ApplyCommand("/process/verbose 1");
// UI->ApplyCommand("/hits/verbose 1");

if (visual) {
  UI->ApplyCommand("/vis/open VRML2FILE");
  // UI->ApplyCommand("/vis/open DAWNFILE");
  UI->ApplyCommand("/vis/scene/create");
  UI->ApplyCommand("/vis/drawVolume");
  if (visuale) {
    UI->ApplyCommand("/vis/scene/add/trajectories");
    //      UI->ApplyCommand("/vis/scene/add/hits");
    UI->ApplyCommand("/vis/sceneHandler/attach");
    UI->ApplyCommand("/tracking/storeTrajectory 1");
  } else {
    UI->ApplyCommand("/vis/viewer/flush");
  }
}
#ifdef G4DISABLE_PARAMETRIZATION
// inactivate the parametrization
UI->ApplyCommand("/process/inactivate G4FastSimulationManagerProcess");
G4cout << "the shower parametrization not used" << G4endl;
#endif

//  UI->ApplyCommand("/geometry/test/grid_test true");

// start a run
// UI->ApplyCommand("/control/suppressAbortion 2");
runManager->SetVerboseLevel(1);
#ifdef G4MYFIT_PARAMETERIZATION
runManager->BeamOn(ParamNumber);
#endif
#ifdef G4MYEM_PARAMETERIZATION
runManager->BeamOn(ParamNumber);
#endif
#if defined(G4MYHA_PARAMETERIZATION) && !defined(G4MYHAMUONS_PARAMETERIZATION)
runManager->BeamOn(ParamNumber);
#endif
#ifndef G4MYFIT_PARAMETERIZATION
#ifndef G4MYEM_PARAMETERIZATION
#if (defined(G4MYHA_PARAMETERIZATION) &&       \
     defined(G4MYHAMUONS_PARAMETERIZATION)) || \
    !defined(G4MYHA_PARAMETERIZATION)
runManager->BeamOn(myGeneratorAction->nevents);
#endif
#endif
#endif

#if (defined(G4MYEM_PARAMETERIZATION) &&    \
     !defined(G4MYK40_PARAMETERIZATION)) || \
    defined(G4MYHA_PARAMETERIZATION)
for (G4int i = 0; i < VertexDistAngleBins; i++) {
  fprintf(outfilePar, "%.20Le %.20Le %.20Le\n", (*myCumNorma)[i],
          (*myCumPhotons)[i], (*myCumPhotonsRms)[i]);
  G4int kstart = i * OMSolidAngleBins;
  G4int kstop = kstart + OMSolidAngleBins;
  for (G4int k = kstart; k < kstop; k++)
    fprintf(outfilePar, "%.20Le %.20Le %.20Le\n", (*myPhotonsTh2)[k],
            (*myPhotonsTh3)[k], (*myPhotonsTh2Th3Num)[k]);
  G4int jstart = i * TimeSolidAngleBins * TimeBins;
  G4int jstop = jstart + TimeSolidAngleBins * TimeBins;
  for (G4int j = jstart; j < jstop; j++)
    fprintf(outfilePar, "%.20Le\n", (*myPhotonsTime)[j]);
}
#endif

if (visual && visuale) {
  UI->ApplyCommand("/vis/viewer/flush");
}
//   session->SessionStart();
// delete session;

// job termination
if (!useANTARESformat)
  fclose(savefile);
else
  delete TheEVTtoWrite;
#ifdef G4MYFIT_PARAMETERIZATION
fclose(outfilePar);
#endif
#if defined(G4MYEM_PARAMETERIZATION) || \
    defined(G4MYHA_PARAMETERIZATION)  // newha
fclose(outfilePar);
#endif

#ifdef G4MYHAMUONS_PARAMETERIZATION
outMuonHAFile->close();
#endif

#ifdef G4VIS_USE
// Visualization, if you choose to have it!
delete visManager;
#endif
delete runManager;
return 0;
}
