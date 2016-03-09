#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>

#include "docopt.h"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4UnitsTable.hh"

#include "KM3Sim.h"
#include "KM3Physics.hh"
#include "KM3PrimaryGeneratorAction.hh"
#include "KM3StackingAction.hh"
#include "KM3TrackingAction.hh"
#include "KM3SteppingAction.hh"
#include "KM3EventAction.hh"
#include "KM3Detector.hh"
/** How to make a simple main:
 *
 * new G4Runmanager
 * rm.SetUnserInit(Detector)
 * rm.SetUnserInit(PhysicsList)
 * rm.SetUnserInit(G4ActionInit)
 * rm.init
 *
 * new uimng
 * uimng.applu(options-verbose)
 *
 * rm.beamon(nevts)
 * del rm
 * return 0
 */

/** how to detector:
 *
 * geometry
 * materials
 * sensitive regions
 * redout schemes of sensitive regions
 */

/** how to physics:
 *
 * particles to be used
 * physics processes tbu
 * range cuts on particles (overrides)
 */

/** how to useraction:
 *
 * add UserActionClasses (see below?)
 * define 1 mandatory UserAction
 */

/* Compilation settings:
 *
 * HADRONIC_COMPILE = True
 * TRACK_INFORMATION = True
 * ENABLE_MIE = True
 * DISABLE_PARAM = True
 * MYFIT_PARAM = False
 * EM_PARAM = False
 * MUON_PARAM = False
 * HA_PARAM = False
 * HAMUON_PARAM = False
 */

static const char USAGE[] =
R"(km3sim.

  Usage:
    km3sim [--seed=<sd>] (-i PARAMS) (-d DETECTOR) (-o OUTFILE)
    km3sim (-h | --help)
    km3sim --version

  Options:
    -i PARAMS         File with physics (seawater etc.) input parameters.
    -d DETECTOR       File with detector geometry.
    -o OUTFILE        File to write events to.
    -h --help         Show this screen.
    --seed=<sd>       Set the RNG seed [default: 42].
    --version         Display the current version.
    --no-mie          Disable mie scattering [default: false]
)";

int main(int argc, char *argv[]) {
  str::map<std::string, doctopt::value> args =
    doctop::docopt(USAGE, {argv + 1, argv + argc}, true, "km3sim v0.1");
  for (auto const &arg : args) {
    std::cout << arg.fist << arg.second << std::endl;
  }
  G4long myseed = atol(argv[2]);
  G4cout << myseed << G4endl;
  CLHEP::HepRandom::setTheSeed(myseed);

  char *Geometry_File = argv[4];
  char *Parameter_File = argv[5];
  G4bool useHEPEvt;
  G4bool useANTARESformat;
  char *fileParticles;
  char *filePythiaParticles;
  FILE *outfilePar;
  G4double ParamEnergy;
  G4int ParamNumber;
  G4int ParamParticle;
  if ((argv[8][0] == 'A') && (argv[8][1] == 'N') && (argv[8][2] == 'T') &&
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
  Mydet->Geometry_File = Geometry_File;
  Mydet->Parameter_File = Parameter_File;
  Mydet->outfilePar = outfilePar;
  runManager->SetUserInitialization(Mydet);

  KM3Physics *MyPhys = new KM3Physics;
  MyPhys->aDetector = Mydet;
  runManager->SetUserInitialization(MyPhys);

  runManager->SetNumberOfEventsToBeStored(0);
  KM3PrimaryGeneratorAction *myGeneratorAction = new KM3PrimaryGeneratorAction;
  myGeneratorAction->fileParticles = fileParticles;
  myGeneratorAction->filePythiaParticles = filePythiaParticles;
  myGeneratorAction->idbeam = ParamParticle;
  myGeneratorAction->outfile = savefile;
  myGeneratorAction->ParamEnergy = ParamEnergy;
  myGeneratorAction->useANTARESformat = useANTARESformat;
  myGeneratorAction->useHEPEvt = useHEPEvt;
  Mydet->MyGenerator = myGeneratorAction;

  KM3TrackingAction *myTracking = new KM3TrackingAction;
  myTracking->TheEVTtoWrite = TheEVTtoWrite;
  myTracking->useANTARESformat = useANTARESformat;
  // link between generator and tracking (to provide number of
  // initial particles to trackingAction
  myGeneratorAction->myTracking = myTracking;
  myGeneratorAction->Initialize();
  runManager->SetUserAction(myGeneratorAction);

  KM3EventAction *event_action = new KM3EventAction;
  event_action->outfile = savefile;
  event_action->TheEVTtoWrite = TheEVTtoWrite;
  event_action->useANTARESformat = useANTARESformat;
  myGeneratorAction->event_action = event_action;
  // generator knows event to set the number of initial particles
  runManager->SetUserAction(event_action);

  Mydet->outfile = savefile;
  Mydet->TheEVTtoWrite = TheEVTtoWrite;
  Mydet->useANTARESformat = useANTARESformat;

  KM3StackingAction *myStacking = new KM3StackingAction;
  KM3SteppingAction *myStepping = new KM3SteppingAction;
  myStacking->SetDetector(Mydet);
  myStepping->myStDetector = Mydet;
  myStepping->event_action = event_action;
  runManager->SetUserAction(myStacking);
  runManager->SetUserAction(myTracking);
  runManager->SetUserAction(myStepping);

  // Initialize G4 kernel
  runManager->Initialize();

  // get the pointer to the UI manager and set verbosities
  G4UImanager *UI = G4UImanager::GetUIpointer();
  G4UIsession *session = 0;
  session = new G4UIterminal();
  // inactivate the parametrization
  UI->ApplyCommand("/process/inactivate G4FastSimulationManagerProcess");

  // start a run
  runManager->SetVerboseLevel(1);
  runManager->BeamOn(myGeneratorAction->nevents);

  // job termination
  if (!useANTARESformat)
    fclose(savefile);
  else
    delete TheEVTtoWrite;

  delete runManager;
  return 0;
}
