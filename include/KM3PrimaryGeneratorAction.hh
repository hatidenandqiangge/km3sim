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
//
//

#ifndef KM3PrimaryGeneratorAction_h
#define KM3PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "KM3TrackingAction.hh"
#include "G4ThreeVector.hh"
#if !defined(G4MYEM_PARAMETERIZATION) && !defined(G4MYHA_PARAMETERIZATION) //newha
#include "KM3EventAction.hh"
#endif

#include <stdio.h>
#include "globals.hh"
#ifdef G4MYMUON_PARAMETERIZATION
#include "KM3MuonParam.hh" //for muon param vs distance
#endif

#ifdef G4MYHAMUONS_PARAMETERIZATION
#include <iostream>
#include <fstream>
#include <iomanip>
#endif

#ifndef G4MYFIT_PARAMETERIZATION
#ifndef G4MYEM_PARAMETERIZATION
#ifndef G4MYHA_PARAMETERIZATION
#ifdef G4HADRONIC_COMPILE

#include "HAVertexMuons.hh"

#endif
#endif
#endif
#endif

class G4Event;
class G4VPrimaryGenerator;

#include "HOURSevtREAD.hh"

class KM3PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:  
  KM3PrimaryGeneratorAction();
  ~KM3PrimaryGeneratorAction();
  
public:
  int nevents;
  FILE* outfile;
  FILE* infile;
  char* filePythiaParticles;
  char* fileParticles;
  G4int numberofParticles;
  void GeneratePrimaries(G4Event* anEvent);
  void Initialize(void);
  G4bool useHEPEvt;
  G4bool useANTARESformat;
  KM3TrackingAction* myTracking;
#if !defined(G4MYEM_PARAMETERIZATION) && !defined(G4MYHA_PARAMETERIZATION) //newha
  KM3EventAction* event_action;
#endif
  G4double ParamEnergy;
  G4int idbeam; //type of injected or produced particles (PDG Code)
  G4double random_R;
  G4ThreeVector position;
  G4ThreeVector direction;
#ifdef G4MYHAMUONS_PARAMETERIZATION
  std::ofstream * outMuonHAFile;
#endif

private:
  G4VPrimaryGenerator* HEPEvt;
  HOURSevtREAD* antaresHEPEvt;
#ifdef G4MYMUON_PARAMETERIZATION
  KM3MuonParam* myMuonParam; //for muon param vs distance
#endif
  G4double EventWeight;
  G4ThreeVector detectorCenter;
  G4double detectorMaxRho;
  G4double detectorMaxz;
  G4double bottomPosition;

#ifndef G4MYFIT_PARAMETERIZATION
#ifndef G4MYEM_PARAMETERIZATION
#ifndef G4MYHA_PARAMETERIZATION
#ifdef G4HADRONIC_COMPILE
  HAVertexMuons * aHAVertexMuons;
#endif
#endif
#endif
#endif

#ifdef G4MYK40_PARAMETERIZATION
  G4double beta(G4double x);
  G4double K40Radius;
#endif

#ifdef G4MYSN_PARAMETERIZATION
  G4double SNRadius;
  G4double NeutrinoTheta,NeutrinoPhi;//momentum vector
#endif

public:
  void PutFromDetector(G4ThreeVector dC, G4double dMR, G4double dMz, G4double bP)
  {detectorCenter=dC;detectorMaxRho=dMR;detectorMaxz=dMz;bottomPosition=bP;};
};

#endif /*KM3PrimaryGeneratorAction_h*/
