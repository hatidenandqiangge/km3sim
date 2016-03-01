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


static const char USAGE[] =
    R"(km3sim.

  Usage:
    km3sim (-h | --help)
    km3sim --version

  Options:
    --seed=<sd>       Set the seed for the random number generator.
    --geometry=<gf>   File with detector geometry.
    -o --output       File to write into.
    -i --input        File with input parameters.
    --em=<ef>         File with ElectoMagnetic parametrization.
    --hadronic=<hf>   File with Hadronic parametrization.
    -h --help         Show this screen.
    --version         Display the current version.
)";

int main(int argc, char *argv[]) {
  str::map<std::string, doctopt::value> args =
    doctop::docopt(USAGE, {argv + 1, argv + argc}, true, "km3sim v0.1");
  for (auto const& arg: args) {
    std::cout << arg.fist << arg.second << std::endl;
  }

  G4long myseed = atol(argv[2]);
  CLHEP::HepRandom::setTheSeed(myseed);

  char *Geometry_File = argv[4];
  char *Parameter_File = argv[5];
  char *EMParametrization_FILE = argv[6];
  char *HAParametrization_FILE = argv[7];

  G4bool useHEPEvt;
  char *fileParticles;
  char *fileParamHAmuons = NULL;
  FILE *outfilePar;
  G4double ParamEnergy;
  G4int ParamNumber;
  G4int ParamParticle;

  //--------------------------------------------------------------------------
  // open the output file
  FILE *savefile;
  HOURSevtWRITE *TheEVTtoWrite;
  TheEVTtoWrite = new HOURSevtWRITE(fileParticles, argv[3]);
  //------------------------------------------------------------------------------------------------------------------------------------------------
  G4RunManager *runManager = new G4RunManager;

  // UserInitialization classes (mandatory)
  KM3Detector *Mydet = new KM3Detector;
  runManager->SetUserInitialization(Mydet);
  KM3Physics *MyPhys = new KM3Physics;
  runManager->SetUserInitialization(MyPhys);
  MyPhys->aDetector = Mydet;
  Mydet->Geometry_File = Geometry_File;
  Mydet->Parameter_File = Parameter_File;
  Mydet->EMParametrization_FILE = EMParametrization_FILE;
  Mydet->HAParametrization_FILE = HAParametrization_FILE;
  Mydet->outfilePar = outfilePar;
#if defined(G4MYEM_PARAMETERIZATION) || defined(G4MYHA_PARAMETERIZATION) // newha
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
  G4int OMSolidAngleBins = 834; // for 834 see KM3SD
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
      for (G4int j = 0; j < TimeBins; j++)
        myPhotonsTime->push_back(zero);
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
  myGeneratorAction->myTracking =
      myTracking; // link between generator and tracking (to provide number of
                  // initial particles to trackingAction
  myGeneratorAction->Initialize();
  runManager->SetUserAction(myGeneratorAction);

  KM3EventAction *event_action = new KM3EventAction;
  runManager->SetUserAction(event_action);
  event_action->outfile = savefile;
  event_action->TheEVTtoWrite = TheEVTtoWrite;
#if !defined(G4MYEM_PARAMETERIZATION) &&                                       \
    !defined(G4MYHA_PARAMETERIZATION) // newha
  myGeneratorAction->event_action = event_action; // generator knows event to
                                                  // set the number of initial
                                                  // particles
#endif
#ifdef G4MYFIT_PARAMETERIZATION
  Mydet->event_action = event_action; // KM3SD must know the momentum of the
                                      // particle when crosses the center of the
                                      // detector
#endif
#if defined(G4MYEM_PARAMETERIZATION) || defined(G4MYHA_PARAMETERIZATION) // newha
  event_action->MyGenerator =
      myGeneratorAction; // event action knows generator for EM param purposes
  event_action->MyStDetector = Mydet;
  event_action->myPhotonsNumber = myPhotonsNumber;
  event_action->myCumPhotons = myCumPhotons;
  event_action->myCumPhotonsRms = myCumPhotonsRms;
  event_action->myCumNorma = myCumNorma;
#endif

  Mydet->outfile = savefile;
  Mydet->TheEVTtoWrite = TheEVTtoWrite;

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

  // get the pointer to the UI manager and set verbosities
  G4UImanager *UI = G4UImanager::GetUIpointer();
  G4UIsession *session = 0;
#ifdef G4UI_USE_TCSH
  session = new G4UIterminal(new G4UItcsh);
#else
  session = new G4UIterminal();
#endif


#ifdef G4DISABLE_PARAMETRIZATION
  // inactivate the parametrization
  UI->ApplyCommand("/process/inactivate G4FastSimulationManagerProcess");
  std::cout << "the shower parametrization not used\n";
#endif

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
#if (defined(G4MYHA_PARAMETERIZATION) &&                                       \
     defined(G4MYHAMUONS_PARAMETERIZATION)) ||                                 \
    !defined(G4MYHA_PARAMETERIZATION)
  runManager->BeamOn(myGeneratorAction->nevents);
#endif
#endif
#endif

#if (defined(G4MYEM_PARAMETERIZATION) &&                                       \
     !defined(G4MYK40_PARAMETERIZATION)) ||                                    \
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

  // job termination
  delete TheEVTtoWrite;
#ifdef G4MYFIT_PARAMETERIZATION
  fclose(outfilePar);
#endif
#if defined(G4MYEM_PARAMETERIZATION) || defined(G4MYHA_PARAMETERIZATION) // newha
  fclose(outfilePar);
#endif

#ifdef G4MYHAMUONS_PARAMETERIZATION
  outMuonHAFile->close();
#endif

  delete runManager;
  return 0;
}
