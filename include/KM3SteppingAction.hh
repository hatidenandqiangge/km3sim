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

#ifndef KM3SteppingAction_h
#define KM3SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "KM3Detector.hh"
#include "KM3EventAction.hh"

class KM3SteppingAction : public G4UserSteppingAction {

public:
  KM3SteppingAction();
  virtual ~KM3SteppingAction(){};
  virtual void UserSteppingAction(const G4Step* aStep);
  KM3Detector* myStDetector;
  KM3EventAction* event_action;
private:
  G4double MuonRange(G4double);
  G4double P7[8];
  G4double P1LOW[2];
  G4double P1HIGH[2];
};

#endif
