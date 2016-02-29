//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
#ifndef KM3StackingAction_H
#define KM3StackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "KM3Detector.hh"
#include "KM3Cherenkov.hh"
#include "KM3EMDeltaFlux.hh"

#ifdef G4MYHAMUONS_PARAMETERIZATION
#include <iostream>
#include <fstream>
#include <iomanip>
#endif

class KM3StackingAction : public G4UserStackingAction
{
public:
  KM3StackingAction();
  virtual ~KM3StackingAction();

public:
  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
  virtual void NewStage();
  virtual void PrepareNewEvent();
  void SetDetector(KM3Detector*);
#ifdef G4MYHAMUONS_PARAMETERIZATION
  std::ofstream * outMuonHAFile;
#endif

private:
  KM3Detector* MyStDetector;
#if !defined(G4MYEM_PARAMETERIZATION) && !defined(G4MYHA_PARAMETERIZATION)
#ifdef G4ENABLE_MIE
#ifndef G4DISABLE_PARAMETRIZATION
  //delta rays parametrization section. Works only for muons primaries (not for showers ve)
  std::vector<G4int> * idprikeep;
  std::vector<G4double> * depenekeep;
  std::vector<G4double> * timekeep;
  std::vector<G4ThreeVector> * poskeep;
  //before  G4int indexkeep;
  KM3EMDeltaFlux* myFlux;
  void CreateAllWaitingPhotons(void);
  ////////////////////////////////////
#endif
#endif
#endif
#ifdef G4MYHAMUONS_PARAMETERIZATION
  void myrotate(G4ThreeVector &x, const G4ThreeVector& p0);
  std::vector <G4double> * theTimes;
  std::vector <G4ThreeVector> * thePositions;
  std::vector <G4ThreeVector> * theMomentums;
#endif

protected :
};


#endif

