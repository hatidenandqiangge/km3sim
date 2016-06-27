#include "KM3EventAction.h"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"

using CLHEP::keV;
using CLHEP::GeV;
using CLHEP::ns;
using CLHEP::m;

void KM3EventAction::BeginOfEventAction(const G4Event *) {
  if (!(G4ParticleTable::GetParticleTable()->GetReadiness())) {
    G4String msg;
    msg = " You are instantiating G4UserEventAction BEFORE your\n";
    msg += "G4VUserPhysicsList is instantiated and assigned to G4RunManager.\n";
    msg +=
        " Such an instantiation is prohibited by Geant4 version 8.0. To fix "
        "this problem,\n";
    msg +=
        "please make sure that your main() instantiates G4VUserPhysicsList "
        "AND\n";
    msg +=
        "set it to G4RunManager before instantiating other user action "
        "classes\n";
    msg += "such as G4UserEventAction.";
    G4Exception("G4UserEventAction::G4UserEventAction()", "Event0001",
                FatalException, msg);
  }
  centerPre.clear();
  centerPost.clear();
  enterPre.clear();
  enterPost.clear();
  leavePre.clear();
  leavePost.clear();
  centerMomentum.clear();
  enterMomentum.clear();
  leaveMomentum.clear();
  centerPosition.clear();
  enterPosition.clear();
  leavePosition.clear();
  centerTime.clear();
  enterTime.clear();
  leaveTime.clear();
  stopPosition.clear();
  stopTime.clear();
#ifdef G4MYMUON_KEEPENERGY
  EnergyAtPosition.clear();
#endif

  G4ThreeVector vzero(0.0, 0.0, 0.0);
  G4double izero;
  izero = 0.0;
  for (G4int ip = 0; ip < numofMuons; ip++) {
    centerPre.push_back(vzero);
    centerPost.push_back(vzero);
    enterPre.push_back(vzero);
    enterPost.push_back(vzero);
    leavePre.push_back(vzero);
    leavePost.push_back(vzero);
    centerMomentum.push_back(izero);
    enterMomentum.push_back(izero);
    leaveMomentum.push_back(izero);
    centerPosition.push_back(vzero);
    enterPosition.push_back(vzero);
    leavePosition.push_back(vzero);
    centerTime.push_back(izero);
    enterTime.push_back(izero);
    leaveTime.push_back(izero);
    stopPosition.push_back(vzero);
    stopTime.push_back(izero);
  }
  TheEVTtoWrite->ReadEvent();
}

void KM3EventAction::EndOfEventAction(const G4Event *) {
  // write the momentums, positions and times to out file
  G4ThreeVector Momentum;
  G4ThreeVector vzero(0.0, 0.0, 0.0);
  G4double izero = 0.0;
  for (G4int ip = 0; ip < numofMuons; ip++) {
    // record entering position//////
    if ((enterPre[ip] != vzero) && (enterPost[ip] != vzero)) {
      Momentum = ((enterMomentum[ip])) * (enterPost[ip] - enterPre[ip]) /
                 (enterPost[ip] - enterPre[ip]).mag();
      Momentum = Momentum.unit();
      TheEVTtoWrite->AddMuonPositionInfo(
          MuonIds[ip], -1, (enterPosition[ip])[0] / m,
          (enterPosition[ip])[1] / m, (enterPosition[ip])[2] / m, Momentum[0],
          Momentum[1], Momentum[2], enterMomentum[ip] / GeV,
          (enterTime[ip]) / ns);
    } else {
      TheEVTtoWrite->AddMuonPositionInfo(MuonIds[ip], -1, 0., 0., 0., 0., 0.,
                                         0., 0., 0.);
    }
    // record center position///////////
    if ((centerPre[ip] != vzero) && (centerPost[ip] != vzero)) {
      Momentum = ((centerMomentum[ip])) * (centerPost[ip] - centerPre[ip]) /
                 (centerPost[ip] - centerPre[ip]).mag();
      Momentum = Momentum.unit();
      TheEVTtoWrite->AddMuonPositionInfo(
          MuonIds[ip], 0, (centerPosition[ip])[0] / m,
          (centerPosition[ip])[1] / m, (centerPosition[ip])[2] / m, Momentum[0],
          Momentum[1], Momentum[2], centerMomentum[ip] / GeV,
          (centerTime[ip]) / ns);
    } else {
      TheEVTtoWrite->AddMuonPositionInfo(MuonIds[ip], 0, 0., 0., 0., 0., 0., 0.,
                                         0., 0.);
    }
    // record leaving position//////////////////////
    if ((leavePre[ip] != vzero) && (leavePost[ip] != vzero)) {
      Momentum = ((leaveMomentum[ip])) * (leavePost[ip] - leavePre[ip]) /
                 (leavePost[ip] - leavePre[ip]).mag();
      Momentum = Momentum.unit();
      TheEVTtoWrite->AddMuonPositionInfo(
          MuonIds[ip], 1, (leavePosition[ip])[0] / m,
          (leavePosition[ip])[1] / m, (leavePosition[ip])[2] / m, Momentum[0],
          Momentum[1], Momentum[2], leaveMomentum[ip] / GeV,
          (leaveTime[ip]) / ns);
    } else {
      TheEVTtoWrite->AddMuonPositionInfo(MuonIds[ip], 1, 0., 0., 0., 0., 0., 0.,
                                         0., 0.);
    }
    // record stopping position//////////////
    TheEVTtoWrite->AddMuonPositionInfo(
        MuonIds[ip], 2, (stopPosition[ip])[0] / m, (stopPosition[ip])[1] / m,
        (stopPosition[ip])[2] / m, (stopTime[ip]) / ns);
  }
// write information of muon energies every 10 meters
#ifdef G4MYMUON_KEEPENERGY
  for (int ien = 0; ien < EnergyAtPosition.size(); ien++)
    EnergyAtPosition[ien] = EnergyAtPosition[ien] / GeV;  // convert to GeV
  TheEVTtoWrite->AddMuonEnergyInfo(EnergyAtPosition);
#endif
  // write to output file
  TheEVTtoWrite->WriteEvent();
}
