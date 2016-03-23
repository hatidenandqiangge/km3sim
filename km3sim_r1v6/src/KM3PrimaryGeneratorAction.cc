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

#include "KM3PrimaryGeneratorAction.hh"

#include "globals.hh"
#include "G4Event.hh"
#include "G4HEPEvtInterface.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

#ifdef G4MYEM_PARAMETERIZATION
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#endif

#ifdef G4MYK40_PARAMETERIZATION
G4double KM3PrimaryGeneratorAction::beta(G4double x) {
  G4double b = sqrt(x * x + 2 * x * 511.) * powf(1300. - x, 2.0) *
               (x + 511.); // see for example hyperphysics
  // "A simple relation for the Fermi function", eq 7
  G4double Z = 20.0;
  G4double E = x + 511.0;
  G4double P = sqrt(E * E - 511.0 * 511.0);
  G4double eta = (Z / 137.0) * E / P;
  G4double fermin =
      fabs(2.0 * 3.14159 * eta / (1.0 - exp(-2.0 * 3.14159 * eta)));
  G4double w = x + 511.0;
  G4double g = Z / 137.0;
  G4double s = sqrt(1.0 - powf(Z / 137.0, 2.0)) - 1.0;
  G4double fermi = fermin * powf(w * w * (1.0 + 4.0 * g * g) - 1.0, s);
  b = b * fermi;
  // half-life determination of K40 by LSC Malonda,Carles, Applied Radiation and
  // Isotopes 56 (2002) 153-156
  w = x / 511. + 1.0;
  G4double wm = 1300. / 511. + 1.0;
  G4double p = sqrt(w * w - 1.0);
  G4double q = wm - w;
  s = powf(p, 6.0) + 7 * powf(p, 4.0) * powf(q, 2.0) +
      7 * powf(p, 2.0) * powf(q, 4.0) + powf(q, 6.0);
  b = b * s;
  return b;
}
#endif

KM3PrimaryGeneratorAction::KM3PrimaryGeneratorAction() {}

KM3PrimaryGeneratorAction::~KM3PrimaryGeneratorAction() {
  if (useHEPEvt && !useANTARESformat)
    delete HEPEvt;
  if (useANTARESformat)
    delete antaresHEPEvt;
#ifdef G4MYMUON_PARAMETERIZATION
  delete myMuonParam; // for muon param vs distance
#endif

#ifndef G4MYFIT_PARAMETERIZATION
#ifndef G4MYEM_PARAMETERIZATION
#ifndef G4MYHA_PARAMETERIZATION
#ifdef G4HADRONIC_COMPILE
#ifndef G4DISABLE_PARAMETRIZATION
  delete aHAVertexMuons;
#endif
#endif
#endif
#endif
#endif
}

void KM3PrimaryGeneratorAction::Initialize() {
  if (useHEPEvt && !useANTARESformat)
    HEPEvt = new G4HEPEvtInterface(filePythiaParticles);
#ifdef G4MYMUON_PARAMETERIZATION
  myMuonParam = new KM3MuonParam(); // for muon param vs distance
#endif
#ifndef G4MYFIT_PARAMETERIZATION
#ifndef G4MYEM_PARAMETERIZATION
#if (defined(G4MYHA_PARAMETERIZATION) &&                                       \
     defined(G4MYHAMUONS_PARAMETERIZATION)) ||                                 \
    !defined(G4MYHA_PARAMETERIZATION)
  if (useANTARESformat) {
    antaresHEPEvt = new HOURSevtREAD(fileParticles);
    nevents = antaresHEPEvt->GetNumberOfEvents();
    useHEPEvt = antaresHEPEvt->IsNeutrinoEvent();
  } else {
    infile = fopen(fileParticles, "r"); // it contains the number of events, the
                                        // particle type, the vertex and
                                        // momentum information
    G4double runtime;
    fscanf(infile, "%d %lf\n", &nevents, &runtime);
  }
#endif
#ifndef G4MYHA_PARAMETERIZATION
#ifdef G4HADRONIC_COMPILE
#ifndef G4DISABLE_PARAMETRIZATION
  aHAVertexMuons =
      new HAVertexMuons("TheMuons_Hadron.dat", "TheMuons_Hadron.index");
#endif
#endif
#endif
#endif
#endif
  if (outfile == NULL && !useANTARESformat)
    G4cout << "ERROR OUTFILE\n" << G4endl;
}
// primary particle generation. Single, multiple (many vertexes) particles and
// neutrino interaction events (single vertex) are supported
// that covers almost everything, except exotic particles (monopoles etc)
void KM3PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {
  static G4int ievent = 0;
  G4int idneu; // type of neutrino interacting (PDG Code)
  G4int idtarget; // type of target if neutrino interaction (PDG Code)
  G4double xneu, yneu, zneu; // neutrino vertex (cm) if neutrino interaction
  G4double pxneu, pyneu,
      pzneu; // neutrino momentum (GeV/c) if neutrino interaction
  G4double xx0, yy0, zz0; // Vertex of injected or produced particles (in cm)
  G4double pxx0, pyy0,
      pzz0; // Momentum of injected or produced particles (in GeV/c)
  G4double t0; // initial time of injected particles(ns)
  ievent++;
#ifndef G4MYFIT_PARAMETERIZATION
#ifndef G4MYEM_PARAMETERIZATION
#if (defined(G4MYHA_PARAMETERIZATION) &&                                       \
     defined(G4MYHAMUONS_PARAMETERIZATION)) ||                                 \
    !defined(G4MYHA_PARAMETERIZATION)
#ifndef G4MYHA_PARAMETERIZATION // newha
  event_action->Initialize();
#endif
  if (ievent == 1 && !useANTARESformat)
    fprintf(outfile, "%d\n", nevents); // write the number of events
  if (!useHEPEvt) {
    idtarget = 0; // the target id is not relevant in case of injected
                  // particles.
    idneu = 0; // neither is the neutrino id
    xneu = 0.0;
    yneu = 0.0;
    zneu = 0.0; // or the neutrino interaction vertex
    pxneu = 0.0;
    pyneu = 0.0;
    pzneu = 0.0; // or the neutrino momentum
    if (!useANTARESformat) {
      fprintf(outfile, "%d %d %d\n", ievent, idneu, idtarget);
      fprintf(outfile, "%.6e %.6e %.6e %.6e %.6e %.6e\n", xneu, yneu, zneu,
              pxneu, pyneu, pzneu);
    }
    if (useANTARESformat) {
      antaresHEPEvt->ReadEvent();
      numberofParticles = antaresHEPEvt->GetNumberOfParticles();
    } else
      fscanf(infile, "%d\n", &numberofParticles);
#ifdef G4MYMUON_PARAMETERIZATION
    // newcode for muon param vs distance // applicable only for simulation of
    // atmospheric muons from EAS
    G4int idbeams[210000];
    G4double t0s[210000], energies[210000], distances[210000];
    std::vector<G4ThreeVector> theMuonsPositions;
    std::vector<G4ThreeVector> theMuonsMomenta;
    std::vector<G4ThreeVector> theMuonsDirections;
    for (G4int ipart = 0; ipart < numberofParticles; ipart++) {
      if (useANTARESformat) {
        antaresHEPEvt->GetParticleInfo(idbeams[ipart], xx0, yy0, zz0, pxx0,
                                       pyy0, pzz0, t0s[ipart]);
      } else
        fscanf(infile, "%d %lf %lf %lf %lf %lf %lf %lf\n", &idbeams[ipart],
               &xx0, &yy0, &zz0, &pxx0, &pyy0, &pzz0, &t0s[ipart]);
      theMuonsPositions.push_back(G4ThreeVector(xx0 * cm, yy0 * cm, zz0 * cm));
      theMuonsMomenta.push_back(
          G4ThreeVector(pxx0 * GeV, pyy0 * GeV, pzz0 * GeV));
      theMuonsDirections.push_back(theMuonsMomenta[ipart].unit());
      t0s[ipart] *= ns;
      energies[ipart] = theMuonsMomenta[ipart].mag();
      // here find the distance to cross the can

      G4ThreeVector distanceV = theMuonsPositions[ipart] - detectorCenter;
      G4double RRR2 = detectorMaxRho * detectorMaxRho;

      distances[ipart] = -1.0;
      // next check if is going to hit the top of the detector
      G4double Ttop = (detectorMaxz - theMuonsPositions[ipart][2]) /
                      theMuonsDirections[ipart][2];
      if (Ttop > 0) {
        G4ThreeVector PointOnTop(
            distanceV[0] + Ttop * theMuonsDirections[ipart][0],
            distanceV[1] + Ttop * theMuonsDirections[ipart][1],
            distanceV[2] + Ttop * theMuonsDirections[ipart][2]);

        G4double dRhoTop =
            PointOnTop[0] * PointOnTop[0] + PointOnTop[1] * PointOnTop[1];
        if (dRhoTop < RRR2)
          distances[ipart] = Ttop;

        // next check if is going to hit the side of the detector
        G4double c =
            distanceV[0] * distanceV[0] + distanceV[1] * distanceV[1] - RRR2;
        if (c > 0.0 && distances[ipart] < 0) {
          G4double a =
              theMuonsDirections[ipart][0] * theMuonsDirections[ipart][0] +
              theMuonsDirections[ipart][1] * theMuonsDirections[ipart][1];
          G4double b = distanceV[0] * theMuonsDirections[ipart][0] +
                       distanceV[1] * theMuonsDirections[ipart][1];
          G4double dia = b * b - a * c;
          if (dia > 0) {
            dia = sqrt(dia);
            G4double SideDist1 = (-b - dia) / a;
            G4double sidez1 = theMuonsPositions[ipart][2] +
                              SideDist1 * theMuonsDirections[ipart][2];
            if ((sidez1 > bottomPosition) && (sidez1 < detectorMaxz) &&
                (SideDist1 > 0))
              distances[ipart] = SideDist1;
          }
        }
      }
    }
    // check if there is any muon that can cross the detector
    G4double distancemax = -1.e9;
    for (G4int ipart = 0; ipart < numberofParticles; ipart++) {
      if (distances[ipart] > distancemax)
        distancemax = distances[ipart];
    }
    if (distancemax <
        0.0) { // no muon that moves to the detector. put a null muon
      numberofParticles = 1;
      EventWeight = 1.0;
      if (!useANTARESformat)
        fprintf(outfile, "%d %.6e\n", numberofParticles, EventWeight);
      G4PrimaryParticle *initialParticle =
          new G4PrimaryParticle(idbeams[0], 0.0, 0.0, 0.01 * GeV);
      G4PrimaryVertex *vertex =
          new G4PrimaryVertex(theMuonsPositions[0], t0s[0]);
      vertex->SetPrimary(initialParticle);
      anEvent->AddPrimaryVertex(vertex);
      G4cout << "Generating vertex with x0= " << theMuonsPositions[0][0]
             << " y0= " << theMuonsPositions[0][1]
             << " z0= " << theMuonsPositions[0][2] << G4endl;
      G4cout << "Generating null particle with px0= 0.0 py0= 0.0 pz0= 10.0MeV "
             << ievent << G4endl;
      // write the beam and target ids (PDG codes), the vertex (cm), momentum
      // (GeV/c)
      if (!useANTARESformat)
        fprintf(outfile, "%d %.6e %.6e %.6e 0.0 0.0 0.01 %.6e\n", idbeams[0],
                theMuonsPositions[0][0] / cm, theMuonsPositions[0][1] / cm,
                theMuonsPositions[0][2] / cm, t0s[0]);
      if ((idbeams[0] == 13) || (idbeams[0] == -13)) {
        event_action->AddPrimaryNumber(1);
      }
    } else { // there are some particles moving to the detector. Consider only
             // the ones that move to.
      for (G4int ipart = 0; ipart < numberofParticles; ipart++) {
        if (distances[ipart] >= 0.0)
          myMuonParam->AddMuon(energies[ipart], distances[ipart]);
      }
      myMuonParam->Initialize();
      // find if there are any capable particles
      G4int numcapable = 0;
      G4int ip = 0;
      for (G4int ipart = 0; ipart < numberofParticles; ipart++) {
        if (distances[ipart] >= 0.0) {
          if (myMuonParam->IsCapable(ip))
            numcapable++;
          ip++;
        }
      }
      if (numcapable ==
          0) { // no muon capable to reach the detector. put a null muon
        numberofParticles = 1;
        EventWeight = 1.0;
        if (!useANTARESformat)
          fprintf(outfile, "%d %.6e\n", numberofParticles, EventWeight);
        G4PrimaryParticle *initialParticle =
            new G4PrimaryParticle(idbeams[0], 0.0, 0.0, 0.01 * GeV);
        G4PrimaryVertex *vertex =
            new G4PrimaryVertex(theMuonsPositions[0], t0s[0]);
        vertex->SetPrimary(initialParticle);
        anEvent->AddPrimaryVertex(vertex);
        G4cout << "Generating vertex with x0= " << theMuonsPositions[0][0]
               << " y0= " << theMuonsPositions[0][1]
               << " z0= " << theMuonsPositions[0][2] << G4endl;
        G4cout
            << "Generating null particle with px0= 0.0 py0= 0.0 pz0= 10.0MeV "
            << ievent << G4endl;
        // write the beam and target ids (PDG codes), the vertex (cm), momentum
        // (GeV/c)
        if (!useANTARESformat)
          fprintf(outfile, "%d %.6e %.6e %.6e 0.0 0.0 0.01 %.6e\n", idbeams[0],
                  theMuonsPositions[0][0] / cm, theMuonsPositions[0][1] / cm,
                  theMuonsPositions[0][2] / cm, t0s[0]);
        if ((idbeams[0] == 13) || (idbeams[0] == -13)) {
          event_action->AddPrimaryNumber(1);
        }
      } else { // there are some muons capable to reach the detector. Consider
               // only the ones that can.
        G4int numberofParticleskeep = numberofParticles;
        do {
          numberofParticles = 0;
          ip = 0;
          for (G4int ipart = 0; ipart < numberofParticleskeep; ipart++) {
            if (distances[ipart] >= 0.0) {
              if (myMuonParam->IsCapable(ip)) {
                distances[ipart] = myMuonParam->GetDistance(ip);
                energies[ipart] = myMuonParam->GetEnergy(ip);
                if (energies[ipart] > 0.0)
                  numberofParticles++;
              } else
                energies[ipart] = -1.0;
              ip++;
            }
          }
        } while (numberofParticles <= 0);
        EventWeight = myMuonParam->GetWeight();
        if (!useANTARESformat)
          fprintf(outfile, "%d %.6e\n", numberofParticles, EventWeight);
        ip = 0;
        for (G4int ipart = 0; ipart < numberofParticleskeep; ipart++) {
          if (distances[ipart] >= 0.0) {
            if (energies[ipart] > 0.0) {
              G4cout << "FERtest " << distances[ipart] << " " << energies[ipart]
                     << G4endl;
              theMuonsMomenta[ipart] =
                  energies[ipart] * theMuonsDirections[ipart];
              theMuonsPositions[ipart] +=
                  distances[ipart] * theMuonsDirections[ipart];
              t0s[ipart] += distances[ipart] / c_light;
              G4PrimaryParticle *initialParticle = new G4PrimaryParticle(
                  idbeams[ipart], theMuonsMomenta[ipart][0],
                  theMuonsMomenta[ipart][1], theMuonsMomenta[ipart][2]);
              G4PrimaryVertex *vertex =
                  new G4PrimaryVertex(theMuonsPositions[ipart], t0s[ipart]);
              vertex->SetPrimary(initialParticle);
              anEvent->AddPrimaryVertex(vertex);
              G4cout << "Generating vertex with x0= "
                     << theMuonsPositions[ipart][0]
                     << " y0= " << theMuonsPositions[ipart][1]
                     << " z0= " << theMuonsPositions[ipart][2] << G4endl;
              G4cout << "Generating particle with px0= "
                     << theMuonsMomenta[ipart][0]
                     << " py0= " << theMuonsMomenta[ipart][1]
                     << " pz0= " << theMuonsMomenta[ipart][2] << "MeV "
                     << ievent << G4endl;
              // write the beam and target ids (PDG codes), the vertex (cm),
              // momentum (GeV/c)
              if (!useANTARESformat)
                fprintf(outfile, "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                        idbeams[ipart], theMuonsPositions[ipart][0] / cm,
                        theMuonsPositions[ipart][1] / cm,
                        theMuonsPositions[ipart][2] / cm,
                        theMuonsMomenta[ipart][0] / GeV,
                        theMuonsMomenta[ipart][1] / GeV,
                        theMuonsMomenta[ipart][2] / GeV, t0s[ipart]);
              if ((idbeams[ipart] == 13) || (idbeams[ipart] == -13)) {
                event_action->AddPrimaryNumber(ip + 1);
              }
              ip++;
            }
          }
        }
      }
      myMuonParam->Finalize(); // stop using the muon distance-energy param
    }
// newcode for muon param vs distance
// ///////////////////////////////////////////////////////////////
#else
    // old code
    // ///////////////////////////////////////////////////////////////////////////////////////
    EventWeight = 1.0;
    if (!useANTARESformat)
      fprintf(outfile, "%d %.6e\n", numberofParticles, EventWeight);
    for (G4int ipart = 0; ipart < numberofParticles; ipart++) {
      // define the particle object and properties from the particle PDG code
      // idbeam
      if (useANTARESformat) {
        antaresHEPEvt->GetParticleInfo(idbeam, xx0, yy0, zz0, pxx0, pyy0, pzz0,
                                       t0);
      } else
        fscanf(infile, "%d %lf %lf %lf %lf %lf %lf %lf\n", &idbeam, &xx0, &yy0,
               &zz0, &pxx0, &pyy0, &pzz0, &t0);
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
      if (!useANTARESformat)
        fprintf(outfile, "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n", idbeam, xx0,
                yy0, zz0, pxx0, pyy0, pzz0, t0);
#ifndef G4MYHA_PARAMETERIZATION // newha
      if ((idbeam == 13) || (idbeam == -13)) {
        event_action->AddPrimaryNumber(ipart + 1);
      }
#endif
    }
// old code
// ///////////////////////////////////////////////////////////////////////////////////////////
#endif
  }
#endif
#endif
#endif
#ifdef G4MYFIT_PARAMETERIZATION
  event_action->Initialize();
  if (!useHEPEvt) {
    numberofParticles = 1;
    idbeam = 13; // muon minus
    G4double startz = -600.0 * meter;
    G4PrimaryParticle *initialParticle =
        new G4PrimaryParticle(idbeam, 0.0, 0.0, ParamEnergy);
    G4PrimaryVertex *vertex =
        new G4PrimaryVertex(G4ThreeVector(0.0, 0.0, startz), 0.0);
    vertex->SetPrimary(initialParticle);
    anEvent->AddPrimaryVertex(vertex);
    event_action->AddPrimaryNumber(1);
    G4cout << "Generating one event " << ievent << G4endl;
  }
#endif
#if defined(G4MYHA_PARAMETERIZATION) && !defined(G4MYHAMUONS_PARAMETERIZATION)
  if (!useHEPEvt) {
    position[0] = 0.0;
    position[1] = 0.0;
    position[2] = 0.0;
    direction[0] = 0.0;
    direction[1] = 0.0;
    direction[2] = 1.0;
    G4double zpos = 3.46; // it is the shift for 100.0GeV pion and kaon zero
                          // long. see also the EM parametrization section in
                          // here
    G4ThreeVector thisPosition = position - zpos * m * direction;
    G4PrimaryVertex *vertex = new G4PrimaryVertex(thisPosition, 0.0);
    numberofParticles = 1;
    // type of incident particle for HA parametrization is input in the command
    // line
    //    idbeam=211; //pion plus
    //    idbeam=321; //kaon plus
    //    idbeam=130; // kaon zero long
    //    idbeam=2112; //Neutron
    //    idbeam=2212; //Proton
    G4PrimaryParticle *initialParticle = new G4PrimaryParticle(
        idbeam, ParamEnergy * direction[0], ParamEnergy * direction[1],
        ParamEnergy * direction[2]);
    vertex->SetPrimary(initialParticle);
    anEvent->AddPrimaryVertex(vertex);
    G4cout << "Generating one event " << ievent << G4endl;
  }
#endif
#ifdef G4MYSN_PARAMETERIZATION
  if (!useHEPEvt) {
    //-----------------antonis-------------------------//
    if (ievent == 1)
      SNRadius = ParamEnergy * meter / GeV;

    static FILE *SNdata = NULL;
    if (SNdata == NULL) {
      SNdata = fopen("data.dat", "r");
      if (SNdata == NULL)
        perror("data.dat");
      fscanf(SNdata, "%lf %lf\n", &NeutrinoTheta, &NeutrinoPhi);
    }
    float NeutrinoEnergy, theta, PositronEnergy;
    fscanf(SNdata, "%f %f %f\n", &NeutrinoEnergy, &PositronEnergy, &theta);
    ParamEnergy = sqrt(powf(PositronEnergy, 2.0) - powf(0.511, 2.0)) * MeV;

    /// direction
    G4double phi = 2.0 * pi * G4UniformRand();
    G4double costheta = cos(theta * deg);
    G4double sintheta = sqrt(1.0 - costheta * costheta);
    G4double cosphi = cos(phi);
    G4double sinphi = sin(phi);
    direction[0] = sintheta * cosphi; // direction x
    direction[1] = sintheta * sinphi; // direction y
    direction[2] = costheta;
    G4ThreeVector NeutrinoDirection;
    if (NeutrinoTheta == -1 && NeutrinoPhi == -1) {
      costheta = 2.0 * G4UniformRand() - 1.0;
      sintheta = sqrt(1.0 - costheta * costheta);
      phi = 2.0 * pi * G4UniformRand();
      cosphi = cos(phi);
      sinphi = sin(phi);
      NeutrinoDirection[0] = sintheta * cosphi;
      NeutrinoDirection[1] = sintheta * sinphi;
      NeutrinoDirection[2] = costheta;
    } else {
      NeutrinoDirection[0] = sin(NeutrinoTheta * deg) * cos(NeutrinoPhi * deg);
      NeutrinoDirection[1] = sin(NeutrinoTheta * deg) * sin(NeutrinoPhi * deg);
      NeutrinoDirection[2] = cos(NeutrinoTheta * deg);
    }

    direction.rotateUz(NeutrinoDirection);

    //////////////POSITION
    G4double rmax;
    do {
      random_R = G4UniformRand() * SNRadius;
      rmax = G4UniformRand() * powf(SNRadius, 2.0);
    } while (random_R * random_R <= rmax || random_R < 21.6 * cm);
    phi = 2.0 * pi * G4UniformRand();
    costheta = 2.0 * G4UniformRand() - 1.0;
    sintheta = sqrt(1.0 - costheta * costheta);
    cosphi = cos(phi);
    sinphi = sin(phi);
    position[0] = random_R * sintheta * cosphi; // position x
    position[1] = random_R * sintheta * sinphi; // position y
    position[2] = random_R * costheta; // position z
    G4PrimaryVertex *vertex = new G4PrimaryVertex(position, 0.0);
    //////////////////////
    numberofParticles = 1;
    G4PrimaryParticle *initialParticle = new G4PrimaryParticle(
        -11,
        ParamEnergy * direction[0], // momentum
        ParamEnergy * direction[1], ParamEnergy * direction[2]);
    vertex->SetPrimary(initialParticle);
    anEvent->AddPrimaryVertex(vertex);
    G4cout << "Generating positron from SN event " << ievent << " with energy "
           << PositronEnergy << G4endl;
    //-----------------antonis-------------------------//
  }
#endif
// laser beam
#ifdef G4MYLASER_PARAMETERIZATION
  if (!useHEPEvt) {
    /// at first initialize the pointers to Rindex,Q_E, glass and gell
    /// transparencies////
    static G4MaterialPropertyVector *QECathod = NULL;
#ifdef G4MY_TRANSPARENCIES
    static G4MaterialPropertyVector *AbsBenth = NULL;
    static G4MaterialPropertyVector *AbsGell = NULL;
#endif
    if (QECathod == NULL) {
      const G4MaterialTable *theMaterialTable = G4Material::GetMaterialTable();
      for (size_t J = 0; J < theMaterialTable->size(); J++) {
        if ((*theMaterialTable)[J]->GetName() == G4String("Cathod")) {
          G4MaterialPropertiesTable *aMaterialPropertiesTable =
              (*theMaterialTable)[J]->GetMaterialPropertiesTable();
          QECathod = aMaterialPropertiesTable->GetProperty("Q_EFF");
        }
#ifdef G4MY_TRANSPARENCIES
        else if ((*theMaterialTable)[J]->GetName() == G4String("Glass")) {
          G4MaterialPropertiesTable *aMaterialPropertiesTable =
              (*theMaterialTable)[J]->GetMaterialPropertiesTable();
          AbsBenth = aMaterialPropertiesTable->GetProperty("ABSLENGTH");
        } else if ((*theMaterialTable)[J]->GetName() == G4String("Gell")) {
          G4MaterialPropertiesTable *aMaterialPropertiesTable =
              (*theMaterialTable)[J]->GetMaterialPropertiesTable();
          AbsGell = aMaterialPropertiesTable->GetProperty("ABSLENGTH");
        }
#endif
      }
    }
    /// end of
    /// initialization//////////////////////////////////////////////////////
    // input: number of photons, wavelength, laser position
    // the DOM is located at (0,0,0).
    G4int Num_Laser_Photons = 100000;
    idbeam = Num_Laser_Photons; // Optical Photon . Here we put the number of
                                // photons for writing purposes on output file
    G4double Wavelength_Laser_Photons = ParamEnergy * nanometer / GeV;
    G4double photonEnergy = h_Planck * c_light / Wavelength_Laser_Photons;
    G4double LaserX, LaserY, LaserZ, LaserDX, LaserDY, LaserDZ;
    if (ievent == 1) {
      FILE *LaserData = fopen("LaserData.dat", "r");
      if (LaserData == NULL)
        perror("LaserData.dat");
      fscanf(LaserData, "%lf %lf %lf %lf %lf %lf\n", &LaserX, &LaserY, &LaserZ,
             &LaserDX, &LaserDY, &LaserDZ);
      fclose(LaserData);
      position[0] = LaserX * m; // laser position x
      position[1] = LaserY * m; // laser position y
      position[2] = LaserZ * m; // laser position z
      direction[0] = LaserDX; // photon direction x
      direction[1] = LaserDY; // photon direction y
      direction[2] = LaserDZ; // photon direction z
      direction = direction.unit(); // normalize the direction vector
      fprintf(outfile, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e\n", LaserX, LaserY,
              LaserZ, LaserDX, LaserDY, LaserDZ, Wavelength_Laser_Photons / nm);
    }

    G4ParticleDefinition *OPDefinition =
        G4OpticalPhoton::OpticalPhotonDefinition();
    G4PrimaryVertex *vertex = new G4PrimaryVertex(position, 0.0);
    numberofParticles = 0;
    for (G4int ip = 0; ip < Num_Laser_Photons;
         ip++) { // Num_Laser_Photons is the number of photns before the
                 // relative QE and transparencies
      G4double qeProb = QECathod->Value(photonEnergy);
#ifdef G4MY_TRANSPARENCIES
      qeProb *= exp(-15.0 / AbsBenth->Value(photonEnergy) -
                    10.0 / AbsGell->Value(photonEnergy));
#endif
      if (G4UniformRand() < qeProb) {
        numberofParticles++;
        G4double px = photonEnergy * direction[0];
        G4double py = photonEnergy * direction[1];
        G4double pz = photonEnergy * direction[2];

        // polarization is normal to direction
        G4double sx = 1.0;
        G4double sy = 0.0;
        G4double sz = 0.0;
        G4ThreeVector photonPolarization(sx, sy, sz);
        photonPolarization.rotateUz(direction);
        sx = photonPolarization[0];
        sy = photonPolarization[1];
        sz = photonPolarization[2];
        //

        // create the photon as primary
        G4PrimaryParticle *initialParticle =
            new G4PrimaryParticle(OPDefinition, px, py, pz);
        initialParticle->SetPolarization(sx, sy, sz);
        vertex->SetPrimary(initialParticle);
      }
      ////////////////////////////////////////////////////////////////////////////
    }
    anEvent->AddPrimaryVertex(vertex);
    G4cout << "Generating one laser pulse " << ievent << " Number of photons "
           << numberofParticles << G4endl;
  }
#endif
// laser beam
#if defined(G4MYK40_PARAMETERIZATION) && !defined(G4MYSN_PARAMETERIZATION) &&  \
    !defined(G4MYLASER_PARAMETERIZATION)
  if (!useHEPEvt) {
    //-----------------antonis-------------------------//
    if (ievent == 1)
      K40Radius = ParamEnergy * meter / GeV;
    position[0] = 0.0; // vertex position x
    position[1] = 0.0; // vertex position y
    position[2] = 0.0; // vertex position z
    G4double phi = 2.0 * pi * G4UniformRand();
    G4double costheta = 2.0 * G4UniformRand() - 1.0;
    G4double sintheta = sqrt(1.0 - costheta * costheta);
    G4double cosphi = cos(phi);
    G4double sinphi = sin(phi);
    direction[0] = sintheta * cosphi; // direction x
    direction[1] = sintheta * sinphi; // direction y
    direction[2] = costheta; // direction z

    G4double rmax;
    do {
      random_R = G4UniformRand() * K40Radius;
      rmax = G4UniformRand() * powf(K40Radius, 2.0);
    } while (random_R * random_R <= rmax || random_R < 21.6 * cm);
    phi = 2.0 * pi * G4UniformRand();
    costheta = 2.0 * G4UniformRand() - 1.0;
    sintheta = sqrt(1.0 - costheta * costheta);
    cosphi = cos(phi);
    sinphi = sin(phi);
    position[0] = random_R * sintheta * cosphi; // position x
    position[1] = random_R * sintheta * sinphi; // position y
    position[2] = random_R * costheta; // position z
    G4PrimaryVertex *vertex = new G4PrimaryVertex(position, 0.0);

    G4double BR = G4UniformRand();
    if (BR <= 0.89338) {
      G4double bmax = 0.24647E15;
      G4double random_T, Tmax;
      do {
        random_T = G4UniformRand() * 1300.0;
        Tmax = G4UniformRand() * bmax;
      } while (beta(random_T) <= Tmax);
      idbeam = 11;
      ParamEnergy = sqrt(random_T * random_T + 2.0 * random_T * 511.0) * keV;
    } else if (0.89338 < BR && BR < 0.89338 + 0.10477) {
      idbeam = 22;
      ParamEnergy = 1460.0 * keV;
    } else {
      idbeam = 22;
      ParamEnergy = 1.0 * keV; // it should be zero, but it core dumps, so put a
                               // small energy (<Ch thres)
    }
    numberofParticles = 1;
    //    idbeam=11; //electron (-11 is positron and 22 is gamma)
    G4PrimaryParticle *initialParticle = new G4PrimaryParticle(
        idbeam,
        ParamEnergy * direction[0], // momentum
        ParamEnergy * direction[1], ParamEnergy * direction[2]);
    vertex->SetPrimary(initialParticle);
    anEvent->AddPrimaryVertex(vertex);
    G4cout << "Generating one K40 event " << ievent << " " << idbeam << " "
           << ParamEnergy << " " << position / m << " " << direction << G4endl;
    //    fprintf(outfile,"%d %.6e %.6e\n",idbeam,ParamEnergy/keV,random_R/m);
    //-----------------antonis-------------------------//
  }
#endif

#if defined(G4MYEM_PARAMETERIZATION) && !defined(G4MYK40_PARAMETERIZATION)
  if (!useHEPEvt) {
    position[0] = 0.0;
    position[1] = 0.0;
    position[2] = 0.0;
    direction[0] = 0.0;
    direction[1] = 0.0;
    direction[2] = 1.0;
    G4PrimaryVertex *vertex;
    if (ParamEnergy > 0.0) {
      // next we tranfer the particle backwards, in order the maximum (median)
      // cherenkov emmision to happen at 0.0
      G4double enepos = log10(ParamEnergy / GeV);
      G4double zpos;
      if (idbeam == 11 || idbeam == -11) {
        if (enepos < 0)
          zpos = 0.90653 + enepos * (0.96676 + enepos * 0.26983);
        else
          zpos = 0.90132 + enepos * 0.81286;
      } else if (idbeam == 22) {
        if (enepos < 0)
          zpos = 1.2245 + enepos * (0.82081 + enepos * 0.20987);
        else
          zpos = 1.2088 + enepos * 0.786;
      }
      zpos *= m;
      G4ThreeVector thisPosition = position - zpos * direction;
      // here we should put the vertex time to -zpos/c_light;
      vertex = new G4PrimaryVertex(thisPosition, 0.0);
    } else {
      vertex = new G4PrimaryVertex(position, 0.0);
    }
    if (ParamEnergy > 0.0) {
      numberofParticles = 1;
      //      idbeam=11; //electron
      G4PrimaryParticle *initialParticle = new G4PrimaryParticle(
          idbeam, ParamEnergy * direction[0], ParamEnergy * direction[1],
          ParamEnergy * direction[2]);
      vertex->SetPrimary(initialParticle);
    } else if (ParamEnergy == -1000.0) { // here is the parametrization for
                                         // delta rays of low energy
      static G4double specin = -0.9832; // this is the (1+a) where a=-1.9832 the
                                        // spectral index of muIoni generated
                                        // electrons (kinetic energy)
      static G4double emin = pow(0.24, specin); // minimum kinetic energy
                                                // 0.24MeV (threshold for
                                                // cherenkov production
      static G4double emax = pow(31.6, specin); // maximum kinetic energy (above
                                                // that em parametrization takes
                                                // over)
      numberofParticles = 10000;
      idbeam = 11; // electron
      for (G4int ip = 0; ip < numberofParticles; ip++) {
        G4double enekin = emin + G4UniformRand() * (emax - emin);
        enekin = pow(enekin, 1.0 / specin);
        G4double SpecMom =
            sqrt((enekin + 0.511) * (enekin + 0.511) -
                 0.511 * 0.511); // 0.511 is the electron mass in MeV
        G4double enelog = log10(enekin);
        G4double costheta =
            0.70351 +
            enelog *
                (0.40702 +
                 enelog *
                     (-0.12039 +
                      enelog *
                          (-0.92576e-1 +
                           enelog * (0.64971e-1 -
                                     0.98130e-2 * enelog)))); // angle vs kinene
        if (costheta > 1.0)
          costheta = 1.0;
        G4double sintheta = sqrt(1.0 - costheta * costheta);
        G4double phi = 2.0 * pi * G4UniformRand();
        G4double cosphi = cos(phi);
        G4double sinphi = sin(phi);
        G4ThreeVector electronMomentum(SpecMom * sintheta * cosphi,
                                       SpecMom * sintheta * sinphi,
                                       SpecMom * costheta);
        electronMomentum.rotateUz(direction);
        G4PrimaryParticle *initialParticle =
            new G4PrimaryParticle(idbeam, electronMomentum[0],
                                  electronMomentum[1], electronMomentum[2]);
        vertex->SetPrimary(initialParticle);
      }
    } else {
      /// at first initialize the pointers to Rindex,Q_E, glass and gell
      /// transparencies////
      static G4MaterialPropertyVector *QECathod = NULL;
      static G4MaterialPropertyVector *Rindex = NULL;
#ifdef G4MY_TRANSPARENCIES
      static G4MaterialPropertyVector *AbsBenth = NULL;
      static G4MaterialPropertyVector *AbsGell = NULL;
#endif
      if (QECathod == NULL) {
        const G4MaterialTable *theMaterialTable =
            G4Material::GetMaterialTable();
        for (size_t J = 0; J < theMaterialTable->size(); J++) {
          if ((*theMaterialTable)[J]->GetName() == G4String("Cathod")) {
            G4MaterialPropertiesTable *aMaterialPropertiesTable =
                (*theMaterialTable)[J]->GetMaterialPropertiesTable();
            QECathod = aMaterialPropertiesTable->GetProperty("Q_EFF");
          } else if ((*theMaterialTable)[J]->GetName() == G4String("Water")) {
            G4MaterialPropertiesTable *aMaterialPropertiesTable =
                (*theMaterialTable)[J]->GetMaterialPropertiesTable();
            Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
          }
#ifdef G4MY_TRANSPARENCIES
          else if ((*theMaterialTable)[J]->GetName() == G4String("Glass")) {
            G4MaterialPropertiesTable *aMaterialPropertiesTable =
                (*theMaterialTable)[J]->GetMaterialPropertiesTable();
            AbsBenth = aMaterialPropertiesTable->GetProperty("ABSLENGTH");
          } else if ((*theMaterialTable)[J]->GetName() == G4String("Gell")) {
            G4MaterialPropertiesTable *aMaterialPropertiesTable =
                (*theMaterialTable)[J]->GetMaterialPropertiesTable();
            AbsGell = aMaterialPropertiesTable->GetProperty("ABSLENGTH");
          }
#endif
        }
      }
      /// end of
      /// initialization//////////////////////////////////////////////////////
      G4ParticleDefinition *OPDefinition =
          G4OpticalPhoton::OpticalPhotonDefinition();
      numberofParticles = 0;
      idbeam = 0; // Optical Photon . is not used , because it has not a
                  // pdgcoding in geant4
      for (G4int ip = 0; ip < 100000;
           ip++) { // 100000 is the number of photns before the QE and
                   // transparencies

        // sample an optical
        // photon//////////////////////////////////////////////
        G4double nMax = Rindex->GetMaxValue();
        G4double maxCos = 1.0 / nMax;
        G4double Pmin = Rindex->GetMinLowEdgeEnergy();
        G4double Pmax = Rindex->GetMaxLowEdgeEnergy();
        G4double dp = Pmax - Pmin;
        G4double maxSin2 = (1.0 - maxCos) * (1.0 + maxCos);

        G4double rand;
        G4double sampledMomentum, sampledRI;
        G4double cosTheta, sin2Theta;

        // sample a phi
        rand = G4UniformRand();
        G4double phi = 2.0 * pi * rand;
        G4double sinPhi = sin(phi);
        G4double cosPhi = cos(phi);

        // Determine photon momentum
        // sample a momentum

        do {
          rand = G4UniformRand();
          sampledMomentum = Pmin + rand * dp;
          sampledRI = Rindex->Value(sampledMomentum);
          cosTheta = 1.0 / sampledRI;

          sin2Theta = (1.0 - cosTheta) * (1.0 + cosTheta);
          rand = G4UniformRand();

        } while (rand * maxSin2 > sin2Theta);
        G4double qeProb = QECathod->Value(sampledMomentum);
#ifdef G4MY_TRANSPARENCIES
        qeProb *= exp(-15.0 / AbsBenth->Value(sampledMomentum) -
                      10.0 / AbsGell->Value(sampledMomentum));
#endif
        if (G4UniformRand() < qeProb) {
          numberofParticles++;
          G4double sinTheta = sqrt(sin2Theta);
          G4double px = sampledMomentum * sinTheta * cosPhi;
          G4double py = sampledMomentum * sinTheta * sinPhi;
          G4double pz = sampledMomentum * cosTheta;

          // Create photon momentum direction vector
          // The momentum direction is still with respect
          // to the coordinate system where the primary
          // particle direction is aligned with the z axis

          G4ThreeVector photonMomentum(px, py, pz);

          // Rotate momentum direction back to global reference
          // system

          photonMomentum.rotateUz(direction);
          // Determine polarization of new photon

          G4double sx = cosTheta * cosPhi;
          G4double sy = cosTheta * sinPhi;
          G4double sz = -sinTheta;

          G4ThreeVector photonPolarization(sx, sy, sz);

          // Rotate back to original coord system

          photonPolarization.rotateUz(direction);
          // create the photon as primary
          G4PrimaryParticle *initialParticle =
              new G4PrimaryParticle(OPDefinition, photonMomentum[0],
                                    photonMomentum[1], photonMomentum[2]);
          initialParticle->SetPolarization(photonPolarization[0],
                                           photonPolarization[1],
                                           photonPolarization[2]);
          vertex->SetPrimary(initialParticle);
        }
        ////////////////////////////////////////////////////////////////////////////
      }
      G4cout << "Number of photons " << numberofParticles << G4endl;
    }
    anEvent->AddPrimaryVertex(vertex);
    G4cout << "Generating one event " << ievent << G4endl;
  }
#endif
#ifndef G4MYFIT_PARAMETERIZATION
#ifndef G4MYEM_PARAMETERIZATION
#ifndef G4MYHA_PARAMETERIZATION
  else {
    t0 = 0.0; // starting particle time is common in neutrino interaction
    if (useANTARESformat) {
      antaresHEPEvt->ReadEvent();
      antaresHEPEvt->GetNeutrinoInfo(idneu, idtarget, xneu, yneu, zneu, pxneu,
                                     pyneu, pzneu, t0);
    } else
      fscanf(infile, "%d %d %lf %lf %lf %lf %lf %lf\n", &idneu, &idtarget,
             &xneu, &yneu, &zneu, &pxneu, &pyneu, &pzneu);
    // Generate the Event (reads from Pythia output file)
    if (useANTARESformat) {
      antaresHEPEvt->GeneratePrimaryVertex(anEvent);
    } else
      HEPEvt->GeneratePrimaryVertex(anEvent);
    // Change the position of the vertex of the event
    anEvent->GetPrimaryVertex(0)->SetPosition(xneu * cm, yneu * cm, zneu * cm);
    G4int numberofVertices = anEvent->GetNumberOfPrimaryVertex();
    numberofParticles = 0;
    for (G4int iv = 0; iv < numberofVertices; iv++)
      numberofParticles += anEvent->GetPrimaryVertex(iv)->GetNumberOfParticle();
    if (!useANTARESformat) {
      fprintf(outfile, "%d %d %d\n", ievent, idneu,
              idtarget); // write the event number (starting with 1) , the beam
                         // and target ids (PDG codes)
      fprintf(outfile, "%.6e %.6e %.6e %.6e %.6e %.6e\n", xneu, yneu, zneu,
              pxneu, pyneu, pzneu); // write the vertex (cm), momentum (GeV/c)
    }
    EventWeight = 1.0;
    if (!useANTARESformat)
      fprintf(outfile, "%d %.6e\n", numberofParticles, EventWeight);
    G4int ipart = 0;
    for (G4int iv = 0; iv < numberofVertices; iv++) {
      for (G4int ip = 0;
           ip < anEvent->GetPrimaryVertex(iv)->GetNumberOfParticle(); ip++) {
        idbeam = anEvent->GetPrimaryVertex(iv)->GetPrimary(ip)->GetPDGcode();
        pxx0 = anEvent->GetPrimaryVertex(iv)->GetPrimary(ip)->GetPx() / GeV;
        pyy0 = anEvent->GetPrimaryVertex(iv)->GetPrimary(ip)->GetPy() / GeV;
        pzz0 = anEvent->GetPrimaryVertex(iv)->GetPrimary(ip)->GetPz() / GeV;
        if (!useANTARESformat)
          fprintf(outfile, "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n", idbeam,
                  xneu, yneu, zneu, pxx0, pyy0, pzz0,
                  t0); // vertex position is common in neutrino interaction
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
#ifdef G4HADRONIC_COMPILE
#ifndef G4DISABLE_PARAMETRIZATION
    // here only if hadronic interactions are used add muons from hadronic
    // vertex
    // but first calculate the energy of hadrons
    direction = G4ThreeVector(0.0, 0.0, 0.0);
    G4double HadronicEnergy = 0.0;
    for (G4int ipart = 0;
         ipart < anEvent->GetPrimaryVertex(0)->GetNumberOfParticle();
         ipart++) { // only the first vertex which is neutrino vertex
      idbeam =
          fabs(anEvent->GetPrimaryVertex(0)->GetPrimary(ipart)->GetPDGcode());
      if (idbeam != 22 && idbeam != 13 && idbeam != 11 && idbeam != 15 &&
          idbeam != 111) { // disregard em particles
        G4ThreeVector thisMomentum =
            anEvent->GetPrimaryVertex(0)->GetPrimary(ipart)->GetMomentum();
        G4double thisMass =
            anEvent->GetPrimaryVertex(0)->GetPrimary(ipart)->GetMass();
        G4double Energy = sqrt(thisMomentum.mag2() + thisMass * thisMass);
        direction += thisMomentum;
        HadronicEnergy += Energy;
      }
    }
    if (HadronicEnergy > 0.0) {
      direction = direction.unit();
      HadronicEnergy -= 938.92 * MeV;
      HadronicEnergy /= GeV;
      G4int NumOfMuonsFromHAVertex =
          aHAVertexMuons->GetNumberOfMuons(HadronicEnergy);
      // create the vertexes and add the muons
      for (G4int iadd = 0; iadd < NumOfMuonsFromHAVertex; iadd++) {
        aHAVertexMuons->ReadMuon();
        // first get the vertex from the class container
        G4ThreeVector aPosition = aHAVertexMuons->GetPosition();
        G4ThreeVector aMomentum = aHAVertexMuons->GetMomentum();
        t0 = aHAVertexMuons->GetTime();
        aPosition.rotateUz(direction);
        aMomentum.rotateUz(direction);
        idbeam = 13; // BUG ??? I have not keep the muon type. Choose in random
        if (G4UniformRand() > 0.5)
          idbeam = -13;
        xx0 = aPosition[0] + xneu;
        yy0 = aPosition[1] + yneu;
        zz0 = aPosition[2] + zneu;
        pxx0 = aMomentum[0];
        pyy0 = aMomentum[1];
        pzz0 = aMomentum[2];
        t0 += anEvent->GetPrimaryVertex(0)->GetT0();
        // define the particle object and properties from the particle PDG code
        // idbeam
        G4PrimaryParticle *initialParticle =
            new G4PrimaryParticle(idbeam, pxx0 * GeV, pyy0 * GeV, pzz0 * GeV);
        G4PrimaryVertex *vertex = new G4PrimaryVertex(
            G4ThreeVector(xx0 * cm, yy0 * cm, zz0 * cm), t0 * ns);
        vertex->SetPrimary(initialParticle);
        anEvent->AddPrimaryVertex(vertex);
      }
      G4cout << "Hadronic run. Event " << ievent << " with hadronic energy "
             << HadronicEnergy << " GeV" << G4endl;
      G4cout << "Number of Initial Particles " << numberofParticles
             << ". Number of aditional muons " << NumOfMuonsFromHAVertex
             << G4endl;
      // add the number of these muons to the numberofParticles
      numberofParticles += NumOfMuonsFromHAVertex;
    }
#endif
#endif
  }
#else
#ifdef G4MYHAMUONS_PARAMETERIZATION
  // Here we do parametrization for hadronic interactions (primary particles
  // that are considered are all except leptons, gamma and pion zero
  // input files must have been generated from dif_cross application with hadron
  // parametrization enabled
  // that means that in the secondaries there are no leptons, gammas and pion
  // zeros
  // the output is just the energetic muons (Ek>1GeV)
  else {
    position[0] = 0.0;
    position[1] = 0.0;
    position[2] = 0.0;
    t0 = 0.0; // starting particle time is common in neutrino interaction
    // Generate the Event (reads from Pythia output file)
    HEPEvt->GeneratePrimaryVertex(anEvent);
    // Change the position of the vertex of the event
    anEvent->GetPrimaryVertex(0)->SetPosition(0.0, 0.0, 0.0);
    numberofParticles = anEvent->GetPrimaryVertex(0)->GetNumberOfParticle();
    direction = G4ThreeVector(0.0, 0.0, 0.0);
    G4float HadronicEnergy = 0.0;
    for (G4int ipart = 0; ipart < numberofParticles; ipart++) {
      idbeam =
          fabs(anEvent->GetPrimaryVertex(0)->GetPrimary(ipart)->GetPDGcode());
      if (idbeam != 22 && idbeam != 13 && idbeam != 11 && idbeam != 15 &&
          idbeam != 111) { // disregard em particles
        G4ThreeVector thisMomentum =
            anEvent->GetPrimaryVertex(0)->GetPrimary(ipart)->GetMomentum();
        G4double thisMass =
            anEvent->GetPrimaryVertex(0)->GetPrimary(ipart)->GetMass();
        G4double Energy = sqrt(thisMomentum.mag2() + thisMass * thisMass);
        direction += thisMomentum;
        HadronicEnergy += (float)Energy;
      }
    }
    direction =
        direction.unit(); // this direction is used in the stacking action
    HadronicEnergy -= 938.92 * MeV;
    HadronicEnergy /= (float)GeV;
    outMuonHAFile->write((char *)&HadronicEnergy, sizeof(HadronicEnergy));
    G4cout << "Hadronic parametrization. Event " << ievent << " with energy "
           << HadronicEnergy << " GeV" << G4endl;
  }
#endif
#endif
#endif
#endif
  myTracking->numofInitialParticles = numberofParticles;
}
