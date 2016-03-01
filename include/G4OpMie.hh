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
// $Id: G4OpMie.hh,v 1.9 2008/07/21 21:08:40 gunter Exp $
// GEANT4 tag $Name: geant4-09-01-patch-01 $
//
//
////////////////////////////////////////////////////////////////////////
// Optical Photon Mie Scattering Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpMie.hh
// Description: Discrete Process -- Mie scattering of optical photons
// Version:     1.0
// Created:     2008-07-21
// Author:      Apostolos Tsirigotis
// Updated:
// mail:        tsirigotis@eap.gr
//
////////////////////////////////////////////////////////////////////////

#ifndef G4OpMie_h
#define G4OpMie_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4OpticalPhoton.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include <vector>
// Class Description:
// Discrete Process -- Mie scattering of optical photons.
// Class inherits publicly from G4VDiscreteProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

struct PhaseFactors {
  double c0;
  double c1;
  double c2;
  double c3;
  double c4;
  double c5;
  double c6;
};

class G4OpMie : public G4VDiscreteProcess {

private:
  //////////////
  // Operators
  //////////////

  // G4OpMie& operator=(const G4OpMie &right);

public: // Without description
        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

  G4OpMie(const G4String &processName = "OpMie", G4ProcessType type = fOptical);

  // G4OpMie(const G4OpMie &right);

  ~G4OpMie();

  ////////////
  // Methods
  ////////////

public: // With description
  G4bool IsApplicable(const G4ParticleDefinition &aParticleType);
  // Returns true -> 'is applicable' only for an optical photon.

  G4double GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition *);
  // Returns the mean free path for Mie scattering in water.
  // --- Not yet implemented for other materials! ---

  G4VParticleChange *PostStepDoIt(const G4Track &aTrack, const G4Step &aStep);
  // This is the method implementing Mie scattering.

private:
  /////////////////////
  // Helper Functions
  /////////////////////
  void BuildThePhysicsTable(void);
  G4double SampleAngle(void);
  G4double PhaseFunction(G4double angle);
  ///////////////////////
  // Class Data Members
  ///////////////////////

private:
  G4int IndexPhaseFunction;
  std::vector<PhaseFactors *> *thePhaseFactors;
  CLHEP::RandGeneral *PhaseRand;
  G4double PhaseArray[1801];
};

////////////////////
// Inline methods
////////////////////

inline G4bool G4OpMie::IsApplicable(const G4ParticleDefinition &aParticleType) {
  return (&aParticleType == G4OpticalPhoton::OpticalPhoton());
}

#endif /* G4OpMie_h */
