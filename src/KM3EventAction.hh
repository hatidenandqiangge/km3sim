//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4UserEventAction.hh,v 1.6 2006/06/29 18:09:17 gunter Exp $
// GEANT4 tag $Name: geant4-08-01-patch-01 $
//
//
//

#ifndef KM3EventAction_h
#define KM3EventAction_h 1
#include "stdio.h"
#include <vector>
#include "G4Types.hh"
#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#if defined(G4MYEM_PARAMETERIZATION) || defined(G4MYHA_PARAMETERIZATION) // newha
#include "KM3PrimaryGeneratorAction.hh"
#include "KM3Detector.hh"
#endif

#include "HOURSevtWRITE.hh"

class G4EventManager;
class G4Event;

// class description:
//
//  This is the base class of one of the user's optional action classes.
// The two methods BeginOfEventAction() and EndOfEventAction() are invoked
// at the beginning and the end of one event processing. These methods are
// invoked by G4EventManager.
//  Be aware that BeginOfEventAction() is invoked when a G4Event object is
// sent to G4EventManager. Thus the primary vertexes/particles have already
// been made by the primary generator. In case the user wants to do something
// before generating primaries (i.e., store random number status), do it in
// the G4VUserPrimaryGeneratorAction concrete class.
//

class KM3EventAction : public G4UserEventAction {
public:
  KM3EventAction() { ; }
  ~KM3EventAction() { ; }
  inline void SetEventManager(G4EventManager *value) { fpEventManager = value; }

public: // with description
  void BeginOfEventAction(const G4Event *anEvent);
  void EndOfEventAction(const G4Event *anEvent);
  // Two virtual method the user can override.
protected:
  G4EventManager *fpEventManager;

public:
  std::vector<G4ThreeVector> centerPre;
  std::vector<G4ThreeVector> centerPost;
  std::vector<G4ThreeVector> enterPre;
  std::vector<G4ThreeVector> enterPost;
  std::vector<G4ThreeVector> leavePre;
  std::vector<G4ThreeVector> leavePost;
  std::vector<double> centerMomentum;
  std::vector<double> enterMomentum;
  std::vector<double> leaveMomentum;
  std::vector<G4ThreeVector> centerPosition;
  std::vector<G4ThreeVector> enterPosition;
  std::vector<G4ThreeVector> leavePosition;
  std::vector<double> centerTime;
  std::vector<double> enterTime;
  std::vector<double> leaveTime;
  std::vector<G4ThreeVector> stopPosition;
  std::vector<double> stopTime;
#ifdef G4MYMUON_KEEPENERGY
  std::vector<double> EnergyAtPosition;
#endif
#if defined(G4MYEM_PARAMETERIZATION) || defined(G4MYHA_PARAMETERIZATION) // newha
  std::vector<long double> *myPhotonsNumber;
  std::vector<long double> *myCumPhotons;
  std::vector<long double> *myCumPhotonsRms;
  std::vector<long double> *myCumNorma;
  KM3PrimaryGeneratorAction *MyGenerator;
  KM3Detector *MyStDetector;
#endif
  bool useANTARESformat;
  HOURSevtWRITE *TheEVTtoWrite;

public:
  FILE *outfile;
  inline void AddPrimaryNumber(int);
  inline int GetSlot(int);
  inline void Initialize(void) { numofMuons = 0; }

private:
  int numofMuons;
  int MuonIds[210000]; // 100 before
};

inline void KM3EventAction::AddPrimaryNumber(int aNumber) {
  MuonIds[numofMuons] = aNumber;
  numofMuons++;
}

inline int KM3EventAction::GetSlot(int aNumber) {
  int slot = -1;
  for (int i = 0; i < numofMuons; i++) {
    if (MuonIds[i] == aNumber) {
      slot = i;
      break;
    }
  }
  return slot;
}

#endif
