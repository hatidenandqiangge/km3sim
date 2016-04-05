#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>

#include "docopt.h"

#include "G4RunManager.h"
#include "G4UImanager.h"
#include "G4UIterminal.h"
#include "G4UItcsh.h"
#include "G4UnitsTable.h"

#include "KM3Sim.h"
#include "KM3Physics.h"
#include "KM3PrimaryGeneratorAction.h"
#include "KM3StackingAction.h"
#include "KM3TrackingAction.h"
#include "KM3SteppingAction.h"
#include "KM3EventAction.h"
#include "KM3Detector.h"
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
      doctop::docopt(USAGE, {argv + 1, argv + argc}, true);

  for (auto const &arg : args) {
    std::cout << arg.fist << arg.second << std::endl;
  }
  G4long myseed = atol(argv[2]);
  G4cout << myseed << G4endl;
  CLHEP::HepRandom::setTheSeed(myseed);

  char *Geometry_File = argv[4];
  char *Parameter_File = argv[5];
  G4bool useHEPEvt;
  G4bool useANTARESformat = true;
  char *fileParticles;
  char *filePythiaParticles;
  G4double ParamEnergy;
  G4int ParamNumber;
  G4int ParamParticle;

  useHEPEvt = true;
  useANTARESformat = true;
  fileParticles = argv[9];

  FILE *savefile;
  HOURSevtWRITE *TheEVTtoWrite;
  TheEVTtoWrite = new HOURSevtWRITE(fileParticles, argv[3]);

  G4RunManager *runManager = new G4RunManager;

  KM3Detector *Mydet = new KM3Detector;
  Mydet->Geometry_File = Geometry_File;
  Mydet->Parameter_File = Parameter_File;
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

  delete TheEVTtoWrite;

  delete runManager;
  return 0;
}
