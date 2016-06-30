#include "KM3PrimaryGeneratorAction.h"

#include "globals.hh"
#include "G4Event.hh"
#include "G4HEPEvtInterface.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

using CLHEP::TeV;
using CLHEP::GeV;
using CLHEP::meter;
using CLHEP::ns;
using CLHEP::cm;
using CLHEP::m;

KM3PrimaryGeneratorAction::KM3PrimaryGeneratorAction() {}

KM3PrimaryGeneratorAction::~KM3PrimaryGeneratorAction() {
  delete antaresHEPEvt;
}

void KM3PrimaryGeneratorAction::Initialize() {
  antaresHEPEvt = new HOURSevtRead(infile_evt);
  nevents = antaresHEPEvt->GetNumberOfEvents();
  useHEPEvt = antaresHEPEvt->IsNeutrinoEvent();
}

// primary particle generation. Single, multiple (many vertexes) particles and
// neutrino interaction events (single vertex) are supported
// that covers almost everything, except exotic particles (monopoles etc)
void KM3PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {
  static G4int ievent = 0;

  // type of neutrino interacting (PDG Code)
  G4int idneu;
  // type of target if neutrino interaction (PDG Code)
  G4int idtarget;
  // neutrino vertex (cm) if neutrino interaction
  G4double xneu, yneu, zneu;
  // neutrino momentum (GeV/c) if neutrino interaction
  G4double pxneu, pyneu, pzneu;
  // Vertex of injected or produced particles (in cm)
  G4double xx0, yy0, zz0;
  // Momentum of injected or produced particles (in GeV/c)
  G4double pxx0, pyy0, pzz0;
  // initial time of injected particles(ns)
  G4double t0;

  ievent++;
  event_action->Initialize();
  if (!useHEPEvt) {
    // the target id is not relevant in case of injected particles.
    idtarget = 0;
    // neither is the neutrino id
    idneu = 0;
    xneu = 0.0;
    yneu = 0.0;
    // or the neutrino interaction vertex
    zneu = 0.0;
    pxneu = 0.0;
    pyneu = 0.0;
    // or the neutrino momentum
    pzneu = 0.0;
    antaresHEPEvt->ReadEvent();
    numberofParticles = antaresHEPEvt->GetNumberOfParticles();

    EventWeight = 1.0;
    for (G4int ipart = 0; ipart < numberofParticles; ipart++) {
      // define the particle object and properties from the particle PDG code
      // idbeam
      antaresHEPEvt->GetParticleInfo(idbeam, xx0, yy0, zz0, pxx0, pyy0, pzz0,
                                     t0);
      G4PrimaryParticle *initialParticle =
          new G4PrimaryParticle(idbeam, pxx0 * GeV, pyy0 * GeV, pzz0 * GeV);
      G4PrimaryVertex *vertex = new G4PrimaryVertex(
          G4ThreeVector(xx0 * cm, yy0 * cm, zz0 * cm), t0 * ns);
      vertex->SetPrimary(initialParticle);
      anEvent->AddPrimaryVertex(vertex);
      G4cout << "Generating vertex with x0=" << xx0 << " y0=" << yy0
             << " z0=" << zz0 << G4endl;
      G4cout << "Generating particle with px0=" << pxx0 << " py0=" << pyy0
             << " pz0=" << pzz0 << " " << ievent << G4endl;
      // write the beam and target ids (PDG codes), the vertex (cm), momentum
      // (GeV/c)
      if ((idbeam == 13) || (idbeam == -13)) {
        event_action->AddPrimaryNumber(ipart + 1);
      }
    }
  } else {
    // starting particle time is common in neutrino interaction
    t0 = 0.0;
    antaresHEPEvt->ReadEvent();
    antaresHEPEvt->GetNeutrinoInfo(idneu, idtarget, xneu, yneu, zneu, pxneu,
                                   pyneu, pzneu, t0);
    // Generate the Event (reads from Pythia output file)
    antaresHEPEvt->GeneratePrimaryVertex(anEvent);
    // Change the position of the vertex of the event
    anEvent->GetPrimaryVertex(0)->SetPosition(xneu * cm, yneu * cm, zneu * cm);
    G4int numberofVertices = anEvent->GetNumberOfPrimaryVertex();
    numberofParticles = 0;
    for (G4int iv = 0; iv < numberofVertices; iv++)
      numberofParticles += anEvent->GetPrimaryVertex(iv)->GetNumberOfParticle();
    EventWeight = 1.0;
    G4int ipart = 0;
    for (G4int iv = 0; iv < numberofVertices; iv++) {
      for (G4int ip = 0;
           ip < anEvent->GetPrimaryVertex(iv)->GetNumberOfParticle(); ip++) {
        idbeam = anEvent->GetPrimaryVertex(iv)->GetPrimary(ip)->GetPDGcode();
        pxx0 = anEvent->GetPrimaryVertex(iv)->GetPrimary(ip)->GetPx() / GeV;
        pyy0 = anEvent->GetPrimaryVertex(iv)->GetPrimary(ip)->GetPy() / GeV;
        pzz0 = anEvent->GetPrimaryVertex(iv)->GetPrimary(ip)->GetPz() / GeV;
        if ((idbeam == 13) || (idbeam == -13)) {
          event_action->AddPrimaryNumber(ipart + 1);
        }
        ipart++;
      }
    }
    G4cout << "Generating one neutrino event: evnum= " << ievent
           << " neutrinotype= " << idneu << " target= " << idtarget
           << " energy= " << sqrt(pxneu * pxneu + pyneu * pyneu + pzneu * pzneu)
           << " GeV" << G4endl;
  }
  myTracking->numofInitialParticles = numberofParticles;
}
