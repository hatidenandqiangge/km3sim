#include <iostream>
#include <fstream>
#include <iomanip>

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
  for (auto const& arg : args) {
    std::cout << arg.fist << arg.second << std::endl;
  }

  G4long myseed = atol(argv[2]);
  CLHEP::HepRandom::setTheSeed(myseed);

  char *Geometry_File = argv[4];
  char *Parameter_File = argv[5];

  bool useHEPEvt;
  char *fileParticles;
  FILE *outfilePar;
  double ParamEnergy;
  int ParamNumber;
  int ParamParticle;

  //--------------------------------------------------------------------------
  // open the output file
  FILE *savefile;
  HOURSevtWRITE *TheEVTtoWrite;
  TheEVTtoWrite = new HOURSevtWRITE(fileParticles, argv[3]);
  //--------------------------------------------------------------------------
  G4RunManager *runManager = new G4RunManager;

  KM3Detector *Mydet = new KM3Detector;
  KM3Physics *MyPhys = new KM3Physics;
  runManager->SetUserInitialization(Mydet);
  runManager->SetUserInitialization(MyPhys);
  MyPhys->aDetector = Mydet;
  Mydet->Geometry_File = Geometry_File;
  Mydet->Parameter_File = Parameter_File;
  Mydet->outfilePar = outfilePar;

  // set mandatory user action class
  runManager->SetNumberOfEventsToBeStored(0);
  // myGeneratorAction and MyTrackingAction
  KM3PrimaryGeneratorAction *myGeneratorAction = new KM3PrimaryGeneratorAction;
  myGeneratorAction->outfile = savefile;
  myGeneratorAction->useHEPEvt = useHEPEvt;
  myGeneratorAction->fileParticles = fileParticles;
  myGeneratorAction->ParamEnergy = ParamEnergy;
  myGeneratorAction->idbeam = ParamParticle;
  Mydet->MyGenerator = myGeneratorAction;

  KM3TrackingAction *myTracking = new KM3TrackingAction;
  myTracking->TheEVTtoWrite = TheEVTtoWrite;
  // link between generator and tracking (to provide number of
  // initial particles to trackingAction
  myGeneratorAction->myTracking = myTracking;
  myGeneratorAction->Initialize();
  runManager->SetUserAction(myGeneratorAction);

  KM3EventAction *event_action = new KM3EventAction;
  runManager->SetUserAction(event_action);
  event_action->outfile = savefile;
  event_action->TheEVTtoWrite = TheEVTtoWrite;
  // generator knows event to set the number of initial particles
  myGeneratorAction->event_action = event_action;

  Mydet->outfile = savefile;
  Mydet->TheEVTtoWrite = TheEVTtoWrite;

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
  std::cout << "the shower parametrization not used\n";

  runManager->SetVerboseLevel(1);
  runManager->BeamOn(ParamNumber);

  // job termination
  delete TheEVTtoWrite;
  delete runManager;
  return 0;
}
