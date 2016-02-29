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
// $Id: ExN05EMShowerModel.hh,v 1.8 2003/06/16 16:49:59 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-patch-01 $
//
// 
//----------------------------------------------
// Parameterisation of e+/e-/gamma producing hits
// The hits are the same as defined in the detailed
// simulation.
//----------------------------------------------
#ifndef KM3HAShowerModel_h
#define KM3HAShowerModel_h 1

#include "math.h"
#include "G4VFastSimulationModel.hh"
#include "G4Step.hh"
#include "G4TouchableHandle.hh"
#include "G4OpticalPhoton.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh" 
#include "G4MaterialPropertiesTable.hh"
#include "KM3Detector.hh"
#include <vector>
#include "KM3HAEnergyFlux.hh"
#include "KM3SD.hh"

class KM3HAShowerModel : public G4VFastSimulationModel
{
public:
  //-------------------------
  // Constructor, destructor
  //-------------------------
  KM3HAShowerModel (G4String, G4Region*);
  KM3HAShowerModel (G4String);
  ~KM3HAShowerModel ();

  //------------------------------
  // Virtual methods of the base
  // class to be coded by the user
  //------------------------------

  // -- IsApplicable
  G4bool IsApplicable(const G4ParticleDefinition&);
  // -- ModelTrigger
  G4bool ModelTrigger(const G4FastTrack &);
  // -- User method DoIt
  void DoIt(const G4FastTrack&, G4FastStep&);
  //  void SetDetector(KM3Detector*);
  KM3Detector* myStDetector;
  KM3SD* aMySD;
  void InitializeFlux(char *,G4double,G4double);

private:
  KM3HAEnergyFlux* myFlux;
  G4double EnergyMin;
  G4double EnergyMax;
  G4double thespeedmaxQE;

};
#endif




