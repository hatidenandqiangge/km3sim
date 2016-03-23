// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02PhysicsList.cc,v 1.8 2000/12/04 16:24:08 maire Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "KM3Physics.hh"

#include "G4UserSpecialCuts.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4ios.hh"
#include "G4FastSimulationManagerProcess.hh" //--apostolis parametrization------------
//#include "g4std/iomanip"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

KM3Physics::KM3Physics() : G4VUserPhysicsList() {
#ifdef G4MYK40_PARAMETERIZATION
  defaultCutValue = 0.001 * mm;
#else
  defaultCutValue = 0.5 * mm;
#endif
  SetVerboseLevel(2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

KM3Physics::~KM3Physics() { delete theCerenkovProcess; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4HADRONIC_COMPILE

#include "G4LeptonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

void KM3Physics::ConstructParticle() {
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program.

  G4BosonConstructor aCBoson;
  G4LeptonConstructor aCLepton;
  G4BaryonConstructor aCBaryon;
  G4MesonConstructor aCMeson;
  G4ShortLivedConstructor aCShort;
  G4IonConstructor aCIon;
  aCBoson.ConstructParticle();
  aCLepton.ConstructParticle();
  aCBaryon.ConstructParticle();
  aCMeson.ConstructParticle();
  aCShort.ConstructParticle();
  aCIon.ConstructParticle();
}

#else

void KM3Physics::ConstructParticle() {
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program.

  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();

  //  e+/-
  G4Electron::ElectronDefinition();
  G4Electron::ElectronDefinition()->SetApplyCutsFlag(true);
  G4Positron::PositronDefinition();
  G4Positron::PositronDefinition()->SetApplyCutsFlag(true);

  // mu+/-
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  // nu_e
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();

  // nu_mu
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  //  mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();

  //  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}

#endif

void KM3Physics::ConstructProcess() {
  AddTransportation();
  AddParameterisation();
  ConstructEM();
#ifdef G4HADRONIC_COMPILE
  ConstructHA(); // construct hadronic processes only in case of Pythia input
                 // (for apparent reasons)
#endif
  ConstructGeneral();
  ConstructOP();

  // tempo list all processes for all particles defined up to the end

  // theParticleIterator->reset();
  // while( (*theParticleIterator)() ){
  //   G4ParticleDefinition* particle = theParticleIterator->value();
  //   G4ProcessManager* pmanager = particle->GetProcessManager();
  //   G4String particleName = particle->GetParticleName();
  //   G4ProcessVector* aVector=pmanager->GetProcessList();
  //   G4int isi=aVector->size();
  //   if(!particle->IsShortLived())
  //     G4cout << "------------------ "<<particleName<<" "<<isi<<"
  //     ------------------"<<G4endl;
  //   else
  //     G4cout << "----SHORTSHORT---- "<<particleName<<" "<<isi<<"
  //     ----SHORTSHORT----"<<G4endl;
  //   for( G4int iii=0 ; iii<isi ; iii++){
  //     G4VProcess* app=(*aVector)[iii];
  //     G4String ass=app->GetProcessName();
  //     G4cout << ass << G4endl;
  //   }
  // }

  // tempo
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
//#include "G4GammaConversionToMuons.hh"
#include "G4PhotoElectricEffect.hh"

// newgeant #include "G4MultipleScattering.hh" //in version 4.9.3.p02 geant4
// threatens than G4MultipleScattering class will be removed in next releases,
// so the next ones are assigned for e+-, mu+-, hadrons and ions, respectively
#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#ifdef G4HADRONIC_COMPILE
#include "G4MuonMinusCaptureAtRest.hh"
#include "G4MuonNuclearProcess.hh" //transition to 4.9.6
#include "G4MuonVDNuclearModel.hh" //transition to 4.9.6
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"
#endif
//#include "G4AnnihiToMuPair.hh"
#include "G4hIonisation.hh"

////#include "G4UserSpecialCuts.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Attention. I must change the high energy of all electromagnetic processes
void KM3Physics::ConstructEM() {
  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition *particle = theParticleIterator->value();
    G4ProcessManager *pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    G4double lowE;
    G4double highE;
    G4int nBins;

    if (particleName == "gamma") {
      // gamma conversion to e+e- pair
      G4GammaConversion *aGammaConversion = new G4GammaConversion();
      lowE = 2 * electron_mass_c2;
      highE = 100.0 * PeV; // 100TeV before
      nBins = 220; // 120 before
      //      aGammaConversion->SetPhysicsTableBining(lowE,highE,nBins);
      //      //not applicable to Geant 4.8
      aGammaConversion->SetLambdaBinning(nBins);
      aGammaConversion->SetMinKinEnergy(lowE);
      aGammaConversion->SetMaxKinEnergy(highE);
      pmanager->AddDiscreteProcess(aGammaConversion);
      // gamma conversion to mu+mu- pair is not likely to happen
      //      pmanager->AddDiscreteProcess(new G4GammaConversionToMuons());
      G4ComptonScattering *aComptonScattering = new G4ComptonScattering();
      lowE = 1 * keV;
      highE = 100.0 * PeV; // 100TeV before
      nBins = 220; // 120 before
      //      aComptonScattering->SetPhysicsTableBining(lowE,highE,nBins);
      //      //not applicable to Geant 4.8
      aComptonScattering->SetLambdaBinning(nBins);
      aComptonScattering->SetMinKinEnergy(lowE);
      aComptonScattering->SetMaxKinEnergy(highE);
      pmanager->AddDiscreteProcess(aComptonScattering);
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
      //      pmanager->AddProcess( new G4UserSpecialCuts(),-1,-1,1);

    } else if (particleName == "e-") {
      // electron
      highE = 100.0 * PeV; // 100TeV before
      nBins = 220; // 120 before
      G4eMultipleScattering *theeminusMultipleScattering =
          new G4eMultipleScattering();
      //      theeminusMultipleScattering->SetBinning(nBins);
      //      theeminusMultipleScattering->SetMaxKinEnergy(highE);
      G4eIonisation *theeminusIonisation = new G4eIonisation();
      theeminusIonisation->SetDEDXBinning(nBins);
      theeminusIonisation->SetLambdaBinning(nBins);
      theeminusIonisation->SetMaxKinEnergy(highE);
      G4eBremsstrahlung *theeminusBremsstrahlung = new G4eBremsstrahlung();
      theeminusBremsstrahlung->SetLambdaBinning(nBins);
      theeminusBremsstrahlung->SetDEDXBinning(nBins);
      theeminusBremsstrahlung->SetMaxKinEnergy(highE);
      //      ((G4eIonisation*)theeminusIonisation)->SetLowerBoundLambda(0.00024);
      //      G4cout<<"lower lambda
      //      e-"<<((G4eIonisation*)theeminusIonisation)->GetLowerBoundLambda()<<G4endl;
      //
      // add processes
      pmanager->AddProcess(theeminusMultipleScattering);
      pmanager->AddProcess(theeminusIonisation);
      pmanager->AddProcess(theeminusBremsstrahlung);
      // pmanager->AddProcess( new G4UserSpecialCuts(),-1,-1,1);
      //
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxAlongStep,
                                   1);
      pmanager->SetProcessOrdering(theeminusIonisation, idxAlongStep, 2);
      // KM3Cherenkov.hh:
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(theeminusIonisation, idxPostStep, 2);
      pmanager->SetProcessOrdering(theeminusBremsstrahlung, idxPostStep, 3);

    } else if (particleName == "e+") {
      // positron
      highE = 100.0 * PeV; // 100TeV before
      nBins = 220; // 120 before
      G4eMultipleScattering *theeplusMultipleScattering =
          new G4eMultipleScattering();
      //      theeplusMultipleScattering->SetBinning(nBins);
      //      theeplusMultipleScattering->SetMaxKinEnergy(highE);
      G4eIonisation *theeplusIonisation = new G4eIonisation();
      theeplusIonisation->SetDEDXBinning(nBins);
      theeplusIonisation->SetLambdaBinning(nBins);
      theeplusIonisation->SetMaxKinEnergy(highE);
      G4eBremsstrahlung *theeplusBremsstrahlung = new G4eBremsstrahlung();
      theeplusBremsstrahlung->SetLambdaBinning(nBins);
      theeplusBremsstrahlung->SetDEDXBinning(nBins);
      theeplusBremsstrahlung->SetMaxKinEnergy(highE);
      G4eplusAnnihilation *theeplusAnnihilation = new G4eplusAnnihilation();
      theeplusAnnihilation->SetLambdaBinning(nBins);
      theeplusAnnihilation->SetMaxKinEnergy(highE);
      //      G4VProcess* theAnnihiToMuPair          = new G4AnnihiToMuPair();
      //      //is not likely to happen
      //      ((G4eIonisation*)theeplusIonisation)->SetLowerBoundLambda(0.00024);
      //      G4cout<<"lower lambda
      //      e+"<<((G4eIonisation*)theeplusIonisation)->GetLowerBoundLambda()<<G4endl;
      //
      // add processes
      pmanager->AddProcess(theeplusMultipleScattering);
      pmanager->AddProcess(theeplusIonisation);

      pmanager->AddProcess(theeplusBremsstrahlung);
      pmanager->AddProcess(theeplusAnnihilation);
      //      pmanager->AddProcess(theAnnihiToMuPair);
      // pmanager->AddProcess( new G4UserSpecialCuts(),-1,-1,1);
      //
      // set ordering for AtRestDoIt
      pmanager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);
      //
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxAlongStep, 1);
      pmanager->SetProcessOrdering(theeplusIonisation, idxAlongStep, 2);
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(theeplusIonisation, idxPostStep, 2);
      pmanager->SetProcessOrdering(theeplusBremsstrahlung, idxPostStep, 3);
      pmanager->SetProcessOrdering(theeplusAnnihilation, idxPostStep, 4);
      //      pmanager->SetProcessOrdering(theAnnihiToMuPair, idxPostStep,5);

    } else if (particleName == "mu+" || particleName == "mu-") {
      // muon
      lowE = 500.0 * GeV;
      highE = 100.0 * PeV; // 100TeV before
      nBins = 220; // 120 before
      G4MuMultipleScattering *aMultipleScattering =
          new G4MuMultipleScattering();
      //      aMultipleScattering->SetBinning(nBins);
      //      aMultipleScattering->SetMaxKinEnergy(highE);
      G4MuBremsstrahlung *aBremsstrahlung = new G4MuBremsstrahlung();
      aBremsstrahlung->SetLambdaBinning(nBins);
      aBremsstrahlung->SetDEDXBinning(nBins);
      aBremsstrahlung->SetMaxKinEnergy(highE);
      G4MuPairProduction *aPairProduction = new G4MuPairProduction();
      aPairProduction->SetLambdaBinning(nBins);
      aPairProduction->SetDEDXBinning(nBins);
      aPairProduction->SetMaxKinEnergy(highE);
      G4MuIonisation *anIonisation = new G4MuIonisation();
      anIonisation->SetLambdaBinning(nBins);
      anIonisation->SetDEDXBinning(nBins);
      anIonisation->SetMaxKinEnergy(highE);

#ifdef G4HADRONIC_COMPILE
      G4MuonNuclearProcess *aMuNuclearInteraction = new G4MuonNuclearProcess();
      aMuNuclearInteraction->SetVerboseLevel(0);
      G4MuonVDNuclearModel *muNucModel = new G4MuonVDNuclearModel();
      muNucModel->SetMaxEnergy(highE);
      aMuNuclearInteraction->RegisterMe(muNucModel);
#endif
      //      ((G4MuPairProduction*)aPairProduction)->SetLowerBoundLambda(0.0021);
      //      G4cout<<"lower lambda pair m
      //      "<<((G4MuPairProduction*)aPairProduction)->GetLowerBoundLambda()<<G4endl;
      //
      // add processes
      pmanager->AddProcess(anIonisation);
      pmanager->AddProcess(aMultipleScattering);
      pmanager->AddProcess(aBremsstrahlung);
      pmanager->AddProcess(aPairProduction);

#ifdef G4HADRONIC_COMPILE
      // new add muon minus capture at rest
      if (particleName == "mu-")
        pmanager->AddRestProcess(new G4MuonMinusCaptureAtRest);
      pmanager->AddDiscreteProcess(aMuNuclearInteraction);
#endif
      //
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep, 1);
      pmanager->SetProcessOrdering(anIonisation, idxAlongStep, 2);
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(anIonisation, idxPostStep, 2);
      pmanager->SetProcessOrdering(aBremsstrahlung, idxPostStep, 3);
      pmanager->SetProcessOrdering(aPairProduction, idxPostStep, 4);
#ifdef G4HADRONIC_COMPILE
      pmanager->SetProcessOrdering(aMuNuclearInteraction, idxPostStep, 5);
#endif

    } else if ((!particle->IsShortLived()) &&
               (particle->GetPDGCharge() != 0.0) &&
               (particle->GetParticleName() != "chargedgeantino")) {
      // all others charged particles except geantino
      highE = 100.0 * PeV; // 100TeV before
      nBins = 220; // 120 before

      G4hMultipleScattering *aMultipleScattering = new G4hMultipleScattering();
      //      aMultipleScattering->SetBinning(nBins);
      //      aMultipleScattering->SetMaxKinEnergy(highE);

      G4hIonisation *anIonisation = new G4hIonisation();
      anIonisation->SetLambdaBinning(nBins);
      anIonisation->SetDEDXBinning(nBins);
      anIonisation->SetMaxKinEnergy(highE);

      ////G4VProcess*  theUserCuts = new G4UserSpecialCuts();

      //
      // add processes
      pmanager->AddProcess(anIonisation);
      pmanager->AddProcess(aMultipleScattering);
      // if( particleName ==  "pi+" ||
      // 	  particleName ==  "pi-" ||
      // 	  particleName ==  "kaon+" ||
      // 	  particleName ==  "kaon-" ||
      // 	  particleName ==  "proton" ||
      // 	  particleName ==  "anti_proton"){
      // 	pmanager->AddProcess(new G4hPairProduction());
      // 	pmanager->AddProcess(new G4hBremsstrahlung());
      // }

      /// pmanager->AddProcess(theUserCuts);

      //
      // set ordering for AtRestDoIt
      ////pmanager->SetProcessOrderingToFirst(theUserCuts,idxAtRest);

      //
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep, 1);
      pmanager->SetProcessOrdering(anIonisation, idxAlongStep, 2);

      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(anIonisation, idxPostStep, 2);

      ////pmanager->SetProcessOrdering(theUserCuts,     idxPostStep,3);
    }
  }
}

#ifdef G4HADRONIC_COMPILE
/////////////////////////////////////////////////////////////////////////////
// Hadronic Physics /////////////////////////////////////////////////////////

// New using physics list QGSP_FTFP_BERT
#include "G4HadronicProcessStore.hh"

#include "G4HadronElasticPhysics.hh"

#include "HadronPhysicsQGSP_FTFP_BERT.hh"

#include "G4AlphaInelasticProcess.hh"
#include "G4LEAlphaInelastic.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4LETritonInelastic.hh"

#include "G4PionMinusAbsorptionAtRest.hh"
#include "G4KaonMinusAbsorption.hh"

void KM3Physics::ConstructHA() {
  G4HadronicProcessStore::Instance()->SetVerbose(0);
  // hadron elastic processes. from G4HadronElasticPhysics
  G4double elimitmaxall = 100.0 * PeV;

  // hadron physics (elastic)
  G4HadronElasticPhysics *theHadronElasticPhysics =
      new G4HadronElasticPhysics(2); // verbose level
  theHadronElasticPhysics->ConstructProcess();

  // hadron physics (inelastic). from HadronPhysicsQGSP_FTFP_BERT builder
  HadronPhysicsQGSP_FTFP_BERT *theQGSP_FTFP_BERT =
      new HadronPhysicsQGSP_FTFP_BERT("hadron", true);
  theQGSP_FTFP_BERT->SetVerboseLevel(0);
  theQGSP_FTFP_BERT->ConstructProcess();

  // put elastic and inelastic models max energy to 100PeV only for barions and
  // mesons.
  // since I dont want to dublicate the QGSP_FTFP_BERT code here and since at
  // this moment
  // I don't have access to HadronicProcessStore, the only way to change the
  // upper limit of the
  // high energy models of many barions and mesons is to access the Process
  // through the particle manager
  // and access the highest energy interaction model of this process through the
  // the G4EnergyRangeManager
  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition *particle = theParticleIterator->value();
    G4String particleName = particle->GetParticleName();
    if (particleName == "anti_lambda" || particleName == "anti_neutron" ||
        particleName == "anti_omega-" || particleName == "anti_proton" ||
        particleName == "anti_sigma+" || particleName == "anti_sigma-" ||
        particleName == "anti_xi-" || particleName == "anti_xi0" ||
        particleName == "kaon+" || particleName == "kaon-" ||
        particleName == "kaon0L" || particleName == "kaon0S" ||
        particleName == "lambda" || particleName == "neutron" ||
        particleName == "omega-" || particleName == "pi+" ||
        particleName == "pi-" || particleName == "proton" ||
        particleName == "sigma+" || particleName == "sigma-" ||
        particleName == "xi-" || particleName == "xi0") {
      G4ProcessManager *pmanager = particle->GetProcessManager();
      G4ProcessVector *aVector = pmanager->GetProcessList();
      G4int isi = aVector->size();
      for (G4int iii = 0; iii < isi; iii++) {
        G4VProcess *app = (*aVector)[iii];
        //	G4cout <<"PPPPPnames "<<particleName<<" "<<
        //app->GetProcessName() <<" "<<app->GetProcessSubType()<<G4endl;
        if (app->GetProcessSubType() == 121 ||
            app->GetProcessSubType() ==
                111) { // is hadronic inelastic or elastic process
          G4HadronicProcess *aHp = (G4HadronicProcess *)(*aVector)[iii];
          G4EnergyRangeManager *aMan = aHp->GetManagerPointer();
          G4cout << aHp->GetProcessName() << G4endl;
          const G4MaterialTable *theMaterialTable =
              G4Material::GetMaterialTable();
          G4Material *aMaterial = NULL;
          const G4Element *anElement = NULL;
          for (size_t J = 0; J < theMaterialTable->size(); J++) {
            if ((*theMaterialTable)[J]->GetName() == G4String("Water")) {
              aMaterial = (*theMaterialTable)[J];
              anElement = aMaterial->GetElement(0);
            }
          }
          if (aMaterial == NULL || anElement == NULL) {
            G4Exception(
                "Material Water and one Element of it not found in ConstructHA",
                "", FatalException, "");
          } else {
            G4HadronicInteraction *anIn =
                aMan->GetHadronicInteraction(99 * TeV, aMaterial, anElement);
            anIn->SetMaxEnergy(elimitmaxall);
            anIn->SetVerboseLevel(0);
          }
        }
      }
    }
  }
  // put inelastic models max energy to 100PeV only for barions and mesons

  // next add also inelastic processes for ions not included in
  // HadronPhysicsQGSP_FTFP_BERT builder
  // only low energy models are available
  // alpha particle
  G4AlphaInelasticProcess *theAlphaInelastic = new G4AlphaInelasticProcess();
  G4ProcessManager *pmanager = G4Alpha::Alpha()->GetProcessManager();
  G4LEAlphaInelastic *theLEAlphaModel = new G4LEAlphaInelastic();
  theAlphaInelastic->RegisterMe(theLEAlphaModel);
  pmanager->AddDiscreteProcess(theAlphaInelastic);
  // Deuteron
  G4DeuteronInelasticProcess *theDeuteronInelastic =
      new G4DeuteronInelasticProcess();
  pmanager = G4Deuteron::Deuteron()->GetProcessManager();
  G4LEDeuteronInelastic *theLEDeuteronModel = new G4LEDeuteronInelastic();
  theDeuteronInelastic->RegisterMe(theLEDeuteronModel);
  pmanager->AddDiscreteProcess(theDeuteronInelastic);
  // Triton
  G4TritonInelasticProcess *theTritonInelastic = new G4TritonInelasticProcess();
  pmanager = G4Triton::Triton()->GetProcessManager();
  G4LETritonInelastic *theLETritonModel = new G4LETritonInelastic();
  theTritonInelastic->RegisterMe(theLETritonModel);
  pmanager->AddDiscreteProcess(theTritonInelastic);

  // next take into account the capture processes for pion and kaon minus
  // PionMinus
  pmanager = G4PionMinus::PionMinus()->GetProcessManager();
  pmanager->AddRestProcess(new G4PionMinusAbsorptionAtRest);

  // KaonMinus
  pmanager = G4KaonMinus::KaonMinus()->GetProcessManager();
  pmanager->AddRestProcess(new G4KaonMinusAbsorption);
}

/* old without using physics
lists--------------------------------------------------------------------

#include "G4ExcitedStringDecay.hh"
#include "G4VLongitudinalStringDecay.hh"
#include "G4BinaryCascade.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSModel.hh"
#include "G4VPartonStringModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4StatMF.hh"
#include "G4FermiBreakUp.hh"
#include "G4Evaporation.hh"
#include "G4TheoFSGenerator.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4HadronElasticProcess.hh"
#include "G4LEXiMinusInelastic.hh"
#include "G4HEXiMinusInelastic.hh"
#include "G4XiMinusInelasticProcess.hh"
#include "G4LEAntiXiMinusInelastic.hh"
#include "G4HEAntiXiMinusInelastic.hh"
#include "G4AntiXiMinusInelasticProcess.hh"
#include "G4LEOmegaMinusInelastic.hh"
#include "G4HEOmegaMinusInelastic.hh"
#include "G4OmegaMinusInelasticProcess.hh"
#include "G4LEAntiOmegaMinusInelastic.hh"
#include "G4HEAntiOmegaMinusInelastic.hh"
#include "G4AntiOmegaMinusInelasticProcess.hh"
#include "G4LEAntiNeutronInelastic.hh"
#include "G4HEAntiNeutronInelastic.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"
#include "G4LELambdaInelastic.hh"
#include "G4HELambdaInelastic.hh"
#include "G4LambdaInelasticProcess.hh"
#include "G4LEAntiLambdaInelastic.hh"
#include "G4HEAntiLambdaInelastic.hh"
#include "G4AntiLambdaInelasticProcess.hh"
#include "G4LESigmaMinusInelastic.hh"
#include "G4HESigmaMinusInelastic.hh"
#include "G4SigmaMinusInelasticProcess.hh"
#include "G4LEAntiSigmaMinusInelastic.hh"
#include "G4HEAntiSigmaMinusInelastic.hh"
#include "G4AntiSigmaMinusInelasticProcess.hh"
#include "G4LESigmaPlusInelastic.hh"
#include "G4HESigmaPlusInelastic.hh"
#include "G4SigmaPlusInelasticProcess.hh"
#include "G4LEAntiSigmaPlusInelastic.hh"
#include "G4HEAntiSigmaPlusInelastic.hh"
#include "G4AntiSigmaPlusInelasticProcess.hh"
#include "G4LEXiZeroInelastic.hh"
#include "G4HEXiZeroInelastic.hh"
#include "G4XiZeroInelasticProcess.hh"
#include "G4LEAntiXiZeroInelastic.hh"
#include "G4HEAntiXiZeroInelastic.hh"
#include "G4AntiXiZeroInelasticProcess.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4KaonMinusAbsorptionAtRest.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4LEAntiProtonInelastic.hh"
#include "G4HEAntiProtonInelastic.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4LCapture.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4LFission.hh"
#include "G4HadronFissionProcess.hh"
#include "G4LElastic.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4PiMinusAbsorptionAtRest.hh"
#include "G4AntiProtonAnnihilationAtRest.hh"
#include "G4QGSMFragmentation.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HEKaonMinusInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4HEProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4HENeutronInelastic.hh"
#include "G4QuasiElasticChannel.hh"
#include "G4ProjectileDiffractiveChannel.hh"

void KM3Physics::ConstructHA()
{
  // ************************************************** //
  // *** preparing inelastic reactions for hadrons *** //
  // ************************************************** //
  //
  // high energy model for proton, neutron, pions and kaons
  G4TheoFSGenerator* theHEModel = new G4TheoFSGenerator; //produces many errors
  // all models for treatment of thermal nucleus
  G4Evaporation* theEvaporation = new G4Evaporation;
  G4FermiBreakUp* theFermiBreakUp = new G4FermiBreakUp;
  G4StatMF* theMF = new G4StatMF;
  // evaporation logic
  G4ExcitationHandler* theHandler = new G4ExcitationHandler;
  theHandler->SetEvaporation(theEvaporation);
  theHandler->SetFermiModel(theFermiBreakUp);
  theHandler->SetMultiFragmentation(theMF);
  theHandler->SetMaxAandZForFermiBreakUp(12, 6);
  theHandler->SetMinEForMultiFrag(3.*MeV);
  // pre-equilibrium stage
  G4PreCompoundModel* thePreEquilib = new G4PreCompoundModel(theHandler);
  thePreEquilib->SetMaxEnergy(70*MeV);
  // a no-cascade generator-precompound interaface
  G4GeneratorPrecompoundInterface* theCascade = new
G4GeneratorPrecompoundInterface;
  theCascade->SetDeExcitation(thePreEquilib);
  // QGSP model
  G4VPartonStringModel* theStringModel = new G4QGSModel<G4QGSParticipants>;
  theHEModel->SetTransport(theCascade);
  theHEModel->SetHighEnergyGenerator(theStringModel);
  //  theHEModel->SetQuasiElasticChannel(new G4QuasiElasticChannel); //include
both but only the first is used
  theHEModel->SetProjectileDiffraction(new G4ProjectileDiffractiveChannel);
  theHEModel->SetMinEnergy(9*TeV); //6GeV before
  theHEModel->SetMaxEnergy(100*PeV); //100TeV before
  // Binary cascade for p,n
  G4BinaryCascade* theCasc = new G4BinaryCascade;
  theCasc->SetMinEnergy(65*MeV);
  theCasc->SetMaxEnergy(6.1*GeV);
  // fragmentation
  G4VLongitudinalStringDecay* theFragmentation = new G4QGSMFragmentation;
  G4ExcitedStringDecay* theStringDecay = new
G4ExcitedStringDecay(theFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);
  //
  // Binary Cascade for Pi
  G4BinaryCascade* theCascForPi = new G4BinaryCascade;
  theCascForPi->SetMinEnergy(0*MeV);
  theCascForPi->SetMaxEnergy(1.5*GeV);


  // *************************** //
  // *** elastic scattering *** //
  // *************************** //
  //
  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  G4LElastic* theElasticModel = new G4LElastic();
  theElasticProcess->RegisterMe(theElasticModel);

  // ***************************************** //
  // *** attaching processes to particles *** //
  // ***************************************** //
  //

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "pi+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      // NOTE: PreCo crahes for Pi+
      // thePionPlusInelasticProcess.RegisterMe(thePreEquilib);
      G4PionPlusInelasticProcess* thePionPlusInelasticProcess = new
G4PionPlusInelasticProcess;
      G4LEPionPlusInelastic* theLEPionPlusInelasticModel = new
G4LEPionPlusInelastic();
      theLEPionPlusInelasticModel->SetMinEnergy(1.4*GeV);
      G4HEPionPlusInelastic* theHEPionPlusInelasticModel = new
G4HEPionPlusInelastic();
      theHEPionPlusInelasticModel->SetMaxEnergy(100*PeV);
      thePionPlusInelasticProcess->RegisterMe(theCascForPi);
      thePionPlusInelasticProcess->RegisterMe(theLEPionPlusInelasticModel);
      thePionPlusInelasticProcess->RegisterMe(theHEPionPlusInelasticModel);
      //      thePionPlusInelasticProcess->RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(thePionPlusInelasticProcess);

    } else if (particleName == "pi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      // thePionMinusInelasticProcess.RegisterMe(thePreEquilib);
      G4PionMinusInelasticProcess* thePionMinusInelasticProcess = new
G4PionMinusInelasticProcess;
      G4LEPionMinusInelastic* theLEPionMinusInelasticModel = new
G4LEPionMinusInelastic();
      theLEPionMinusInelasticModel->SetMinEnergy(1.4*GeV);
      G4HEPionMinusInelastic* theHEPionMinusInelasticModel = new
G4HEPionMinusInelastic();
      theHEPionMinusInelasticModel->SetMaxEnergy(100*PeV);
      thePionMinusInelasticProcess->RegisterMe(theCascForPi);
      thePionMinusInelasticProcess->RegisterMe(theLEPionMinusInelasticModel);
      thePionMinusInelasticProcess->RegisterMe(theHEPionMinusInelasticModel);
      //      thePionMinusInelasticProcess->RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(thePionMinusInelasticProcess);
      pmanager->AddRestProcess(new G4PiMinusAbsorptionAtRest, ordDefault);

    } else if (particleName == "kaon+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LEKaonPlusInelastic* theLEKaonPlusInelasticModel = new
G4LEKaonPlusInelastic();
      G4HEKaonPlusInelastic* theHEKaonPlusInelasticModel = new
G4HEKaonPlusInelastic();
      theHEKaonPlusInelasticModel->SetMaxEnergy(100*PeV);
      G4KaonPlusInelasticProcess* theKaonPlusInelasticProcess = new
G4KaonPlusInelasticProcess;
      theKaonPlusInelasticProcess->RegisterMe(theLEKaonPlusInelasticModel);
      theKaonPlusInelasticProcess->RegisterMe(theHEKaonPlusInelasticModel);
      //      theKaonPlusInelasticProcess->RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(theKaonPlusInelasticProcess);

    } else if (particleName == "kaon0S") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LEKaonZeroSInelastic* theLEKaonZeroSInelasticModel = new
G4LEKaonZeroSInelastic();
      G4HEKaonZeroInelastic* theHEKaonZeroInelasticModel = new
G4HEKaonZeroInelastic();
      theHEKaonZeroInelasticModel->SetMaxEnergy(100*PeV);
      G4KaonZeroSInelasticProcess* theKaonZeroSInelasticProcess =new
G4KaonZeroSInelasticProcess;
      theKaonZeroSInelasticProcess->RegisterMe(theLEKaonZeroSInelasticModel);
      theKaonZeroSInelasticProcess->RegisterMe(theHEKaonZeroInelasticModel);
      //      theKaonZeroSInelasticProcess->RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(theKaonZeroSInelasticProcess);

    } else if (particleName == "kaon0L") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LEKaonZeroLInelastic* theLEKaonZeroLInelasticModel = new
G4LEKaonZeroLInelastic();
      G4HEKaonZeroInelastic* theHEKaonZeroInelasticModel = new
G4HEKaonZeroInelastic();
      theHEKaonZeroInelasticModel->SetMaxEnergy(100*PeV);
      G4KaonZeroLInelasticProcess* theKaonZeroLInelasticProcess = new
G4KaonZeroLInelasticProcess;
      theKaonZeroLInelasticProcess->RegisterMe(theLEKaonZeroLInelasticModel);
      theKaonZeroLInelasticProcess->RegisterMe(theHEKaonZeroInelasticModel);
      //      theKaonZeroLInelasticProcess->RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(theKaonZeroLInelasticProcess);

    } else if (particleName == "kaon-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LEKaonMinusInelastic* theLEKaonMinusInelasticModel = new
G4LEKaonMinusInelastic();
      G4HEKaonMinusInelastic* theHEKaonMinusInelasticModel = new
G4HEKaonMinusInelastic();
      theHEKaonMinusInelasticModel->SetMaxEnergy(100*PeV);
      G4KaonMinusInelasticProcess* theKaonMinusInelasticProcess = new
G4KaonMinusInelasticProcess;
      theKaonMinusInelasticProcess->RegisterMe(theLEKaonMinusInelasticModel);
      theKaonMinusInelasticProcess->RegisterMe(theHEKaonMinusInelasticModel);
      //      theKaonMinusInelasticProcess->RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(theKaonMinusInelasticProcess);
      pmanager->AddRestProcess(new G4KaonMinusAbsorptionAtRest, ordDefault);

    } else if (particleName == "proton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LEProtonInelastic* theLEProtonInelasticModel = new
G4LEProtonInelastic();
      theLEProtonInelasticModel->SetMinEnergy(6.0*GeV);
      G4HEProtonInelastic* theHEProtonInelasticModel = new
G4HEProtonInelastic();
      theHEProtonInelasticModel->SetMaxEnergy(100*PeV);
      G4ProtonInelasticProcess* theProtonInelasticProcess = new
G4ProtonInelasticProcess;
      theProtonInelasticProcess->RegisterMe(thePreEquilib);
      theProtonInelasticProcess->RegisterMe(theCasc);
      theProtonInelasticProcess->RegisterMe(theLEProtonInelasticModel);
      theProtonInelasticProcess->RegisterMe(theHEProtonInelasticModel);
      //      theProtonInelasticProcess->RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(theProtonInelasticProcess);

    } else if (particleName == "anti_proton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LEAntiProtonInelastic* theLEAntiProtonInelasticModel = new
G4LEAntiProtonInelastic();
      G4HEAntiProtonInelastic* theHEAntiProtonInelasticModel = new
G4HEAntiProtonInelastic();
      theHEAntiProtonInelasticModel->SetMaxEnergy(100*PeV);
      G4AntiProtonInelasticProcess* theAntiProtonInelasticProcess = new
G4AntiProtonInelasticProcess;
      theAntiProtonInelasticProcess->RegisterMe(theLEAntiProtonInelasticModel);
      theAntiProtonInelasticProcess->RegisterMe(theHEAntiProtonInelasticModel);
      //      theAntiProtonInelasticProcess->RegisterMe(theHEModel); //it was
missing before
      pmanager->AddDiscreteProcess(theAntiProtonInelasticProcess);
      pmanager->AddRestProcess( new G4AntiProtonAnnihilationAtRest, ordDefault);

    } else if (particleName == "neutron") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LENeutronInelastic* theLENeutronInelasticModel = new
G4LENeutronInelastic();
      theLENeutronInelasticModel->SetMinEnergy(6.0*GeV);
      G4HENeutronInelastic* theHENeutronInelasticModel = new
G4HENeutronInelastic();
      theHENeutronInelasticModel->SetMaxEnergy(100*PeV);
      G4NeutronInelasticProcess* theNeutronInelasticProcess = new
G4NeutronInelasticProcess;
      theNeutronInelasticProcess->RegisterMe(thePreEquilib);
      theNeutronInelasticProcess->RegisterMe(theCasc);
      theNeutronInelasticProcess->RegisterMe(theLENeutronInelasticModel);
      theNeutronInelasticProcess->RegisterMe(theHENeutronInelasticModel);
      //      theNeutronInelasticProcess->RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(theNeutronInelasticProcess);
      // capture
      G4LCapture* theNeutronCaptureModel1 = new G4LCapture();
      //   theNeutronCaptureModel1->SetMinEnergy(19*MeV);
      G4HadronCaptureProcess* theNeutronCaptureProcess = new
G4HadronCaptureProcess;
      theNeutronCaptureProcess->RegisterMe(theNeutronCaptureModel1);
      //   theNeutronCaptureModel2 = new G4NeutronHPCapture;
      //   theNeutronCaptureProcess.RegisterMe(theNeutronCaptureModel2);
      //   theNeutronCaptureData = new G4NeutronHPCaptureData;
      //   theNeutronCaptureProcess.AddDataSet(theNeutronCaptureData);
      pmanager->AddDiscreteProcess(theNeutronCaptureProcess);
      // fission
      G4LFission* theNeutronFissionModel = new G4LFission();
      G4HadronFissionProcess* theNeutronFissionProcess = new
G4HadronFissionProcess;
      theNeutronFissionProcess->RegisterMe(theNeutronFissionModel);
      pmanager->AddDiscreteProcess(theNeutronFissionProcess);

    } else if (particleName == "anti_neutron") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LEAntiNeutronInelastic* theLEAntiNeutronInelasticModel = new
G4LEAntiNeutronInelastic();
      G4HEAntiNeutronInelastic* theHEAntiNeutronInelasticModel = new
G4HEAntiNeutronInelastic();
      theHEAntiNeutronInelasticModel->SetMaxEnergy(100*PeV);
      G4AntiNeutronInelasticProcess* theAntiNeutronInelasticProcess = new
G4AntiNeutronInelasticProcess;
      theAntiNeutronInelasticProcess->RegisterMe(theLEAntiNeutronInelasticModel);
      theAntiNeutronInelasticProcess->RegisterMe(theHEAntiNeutronInelasticModel);
      //      theAntiNeutronInelasticProcess->RegisterMe(theHEModel); //it was
missing before
      pmanager->AddDiscreteProcess(theAntiNeutronInelasticProcess);
      pmanager->AddRestProcess(new G4AntiNeutronAnnihilationAtRest,ordDefault);

    } else if (particleName == "lambda") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LELambdaInelastic* theLELambdaInelasticModel = new
G4LELambdaInelastic();
      G4HELambdaInelastic* theHELambdaInelasticModel = new
G4HELambdaInelastic();
      theHELambdaInelasticModel->SetMaxEnergy(100*PeV);
      G4LambdaInelasticProcess* theLambdaInelasticProcess = new
G4LambdaInelasticProcess;
      theLambdaInelasticProcess->RegisterMe(theLELambdaInelasticModel);
      theLambdaInelasticProcess->RegisterMe(theHELambdaInelasticModel);
      //      theLambdaInelasticProcess->RegisterMe(theHEModel); //it was
missing before
      pmanager->AddDiscreteProcess(theLambdaInelasticProcess);

    } else if (particleName == "anti_lambda") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LEAntiLambdaInelastic* theLEAntiLambdaInelasticModel = new
G4LEAntiLambdaInelastic();
      G4HEAntiLambdaInelastic* theHEAntiLambdaInelasticModel = new
G4HEAntiLambdaInelastic();
      theHEAntiLambdaInelasticModel->SetMaxEnergy(100*PeV);
      G4AntiLambdaInelasticProcess* theAntiLambdaInelasticProcess = new
G4AntiLambdaInelasticProcess;
      theAntiLambdaInelasticProcess->RegisterMe(theLEAntiLambdaInelasticModel);
      theAntiLambdaInelasticProcess->RegisterMe(theHEAntiLambdaInelasticModel);
      //      theAntiLambdaInelasticProcess->RegisterMe(theHEModel); //it was
missing before
      pmanager->AddDiscreteProcess(theAntiLambdaInelasticProcess);

    } else if (particleName == "sigma-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LESigmaMinusInelastic* theLESigmaMinusInelasticModel = new
G4LESigmaMinusInelastic();
      G4HESigmaMinusInelastic* theHESigmaMinusInelasticModel = new
G4HESigmaMinusInelastic();
      theHESigmaMinusInelasticModel->SetMaxEnergy(100*PeV);
      G4SigmaMinusInelasticProcess* theSigmaMinusInelasticProcess = new
G4SigmaMinusInelasticProcess;
      theSigmaMinusInelasticProcess->RegisterMe(theLESigmaMinusInelasticModel);
      theSigmaMinusInelasticProcess->RegisterMe(theHESigmaMinusInelasticModel);
      //      theSigmaMinusInelasticProcess->RegisterMe(theHEModel); //it was
missing before
      pmanager->AddDiscreteProcess(theSigmaMinusInelasticProcess);

    } else if (particleName == "anti_sigma-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LEAntiSigmaMinusInelastic* theLEAntiSigmaMinusInelasticModel = new
G4LEAntiSigmaMinusInelastic();
      G4HEAntiSigmaMinusInelastic* theHEAntiSigmaMinusInelasticModel = new
G4HEAntiSigmaMinusInelastic();
      theHEAntiSigmaMinusInelasticModel->SetMaxEnergy(100*PeV);
      G4AntiSigmaMinusInelasticProcess* theAntiSigmaMinusInelasticProcess = new
G4AntiSigmaMinusInelasticProcess;
      theAntiSigmaMinusInelasticProcess->RegisterMe(theLEAntiSigmaMinusInelasticModel);
      theAntiSigmaMinusInelasticProcess->RegisterMe(theHEAntiSigmaMinusInelasticModel);
      //      theAntiSigmaMinusInelasticProcess->RegisterMe(theHEModel); //it
was missing before
      pmanager->AddDiscreteProcess(theAntiSigmaMinusInelasticProcess);

    } else if (particleName == "sigma+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LESigmaPlusInelastic* theLESigmaPlusInelasticModel = new
G4LESigmaPlusInelastic();
      G4HESigmaPlusInelastic* theHESigmaPlusInelasticModel = new
G4HESigmaPlusInelastic();
      theHESigmaPlusInelasticModel->SetMaxEnergy(100*PeV);
      G4SigmaPlusInelasticProcess* theSigmaPlusInelasticProcess = new
G4SigmaPlusInelasticProcess;
      theSigmaPlusInelasticProcess->RegisterMe(theLESigmaPlusInelasticModel);
      theSigmaPlusInelasticProcess->RegisterMe(theHESigmaPlusInelasticModel);
      //      theSigmaPlusInelasticProcess->RegisterMe(theHEModel); //it was
missing before
      pmanager->AddDiscreteProcess(theSigmaPlusInelasticProcess);

    } else if (particleName == "anti_sigma+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LEAntiSigmaPlusInelastic* theLEAntiSigmaPlusInelasticModel = new
G4LEAntiSigmaPlusInelastic();
      G4HEAntiSigmaPlusInelastic* theHEAntiSigmaPlusInelasticModel = new
G4HEAntiSigmaPlusInelastic();
      theHEAntiSigmaPlusInelasticModel->SetMaxEnergy(100*PeV);
      G4AntiSigmaPlusInelasticProcess* theAntiSigmaPlusInelasticProcess = new
G4AntiSigmaPlusInelasticProcess;
      theAntiSigmaPlusInelasticProcess->RegisterMe(theLEAntiSigmaPlusInelasticModel);
      theAntiSigmaPlusInelasticProcess->RegisterMe(theHEAntiSigmaPlusInelasticModel);
      //      theAntiSigmaPlusInelasticProcess->RegisterMe(theHEModel); //it was
missing before
      pmanager->AddDiscreteProcess(theAntiSigmaPlusInelasticProcess);

    } else if (particleName == "xi0") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LEXiZeroInelastic* theLEXiZeroInelasticModel = new
G4LEXiZeroInelastic();
      G4HEXiZeroInelastic* theHEXiZeroInelasticModel = new
G4HEXiZeroInelastic();
      theHEXiZeroInelasticModel->SetMaxEnergy(100*PeV);
      G4XiZeroInelasticProcess* theXiZeroInelasticProcess = new
G4XiZeroInelasticProcess;
      theXiZeroInelasticProcess->RegisterMe(theLEXiZeroInelasticModel);
      theXiZeroInelasticProcess->RegisterMe(theHEXiZeroInelasticModel);
      //      theXiZeroInelasticProcess->RegisterMe(theHEModel); //it was
missing before
      pmanager->AddDiscreteProcess(theXiZeroInelasticProcess);

    } else if (particleName == "anti_xi0") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LEAntiXiZeroInelastic* theLEAntiXiZeroInelasticModel = new
G4LEAntiXiZeroInelastic();
      G4HEAntiXiZeroInelastic* theHEAntiXiZeroInelasticModel = new
G4HEAntiXiZeroInelastic();
      theHEAntiXiZeroInelasticModel->SetMaxEnergy(100*PeV);
      G4AntiXiZeroInelasticProcess* theAntiXiZeroInelasticProcess = new
G4AntiXiZeroInelasticProcess;
      theAntiXiZeroInelasticProcess->RegisterMe(theLEAntiXiZeroInelasticModel);
      theAntiXiZeroInelasticProcess->RegisterMe(theHEAntiXiZeroInelasticModel);
      //      theAntiXiZeroInelasticProcess->RegisterMe(theHEModel); //it was
missing before
      pmanager->AddDiscreteProcess(theAntiXiZeroInelasticProcess);

    } else if (particleName == "xi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LEXiMinusInelastic* theLEXiMinusInelasticModel = new
G4LEXiMinusInelastic();
      G4HEXiMinusInelastic* theHEXiMinusInelasticModel = new
G4HEXiMinusInelastic();
      theHEXiMinusInelasticModel->SetMaxEnergy(100*PeV);
      G4XiMinusInelasticProcess* theXiMinusInelasticProcess = new
G4XiMinusInelasticProcess;
      theXiMinusInelasticProcess->RegisterMe(theLEXiMinusInelasticModel);
      theXiMinusInelasticProcess->RegisterMe(theHEXiMinusInelasticModel);
      //      theXiMinusInelasticProcess->RegisterMe(theHEModel); //it was
missing before
      pmanager->AddDiscreteProcess(theXiMinusInelasticProcess);

    } else if (particleName == "anti_xi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LEAntiXiMinusInelastic* theLEAntiXiMinusInelasticModel = new
G4LEAntiXiMinusInelastic();
      G4HEAntiXiMinusInelastic* theHEAntiXiMinusInelasticModel = new
G4HEAntiXiMinusInelastic();
      theHEAntiXiMinusInelasticModel->SetMaxEnergy(100*PeV);
      G4AntiXiMinusInelasticProcess* theAntiXiMinusInelasticProcess = new
G4AntiXiMinusInelasticProcess;
      theAntiXiMinusInelasticProcess->RegisterMe(theLEAntiXiMinusInelasticModel);
      theAntiXiMinusInelasticProcess->RegisterMe(theHEAntiXiMinusInelasticModel);
      //      theAntiXiMinusInelasticProcess->RegisterMe(theHEModel); //it was
missing before
      pmanager->AddDiscreteProcess(theAntiXiMinusInelasticProcess);
    } else if (particleName == "omega-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LEOmegaMinusInelastic* theLEOmegaMinusInelasticModel = new
G4LEOmegaMinusInelastic();
      G4HEOmegaMinusInelastic* theHEOmegaMinusInelasticModel = new
G4HEOmegaMinusInelastic();
      theHEOmegaMinusInelasticModel->SetMaxEnergy(100*PeV);
      G4OmegaMinusInelasticProcess* theOmegaMinusInelasticProcess = new
G4OmegaMinusInelasticProcess;
      theOmegaMinusInelasticProcess->RegisterMe(theLEOmegaMinusInelasticModel);
      theOmegaMinusInelasticProcess->RegisterMe(theHEOmegaMinusInelasticModel);
      //      theOmegaMinusInelasticProcess->RegisterMe(theHEModel); //it was
missing before
      pmanager->AddDiscreteProcess(theOmegaMinusInelasticProcess);

    } else if (particleName == "anti_omega-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LEAntiOmegaMinusInelastic* theLEAntiOmegaMinusInelasticModel = new
G4LEAntiOmegaMinusInelastic();
      G4HEAntiOmegaMinusInelastic* theHEAntiOmegaMinusInelasticModel = new
G4HEAntiOmegaMinusInelastic();
      theHEAntiOmegaMinusInelasticModel->SetMaxEnergy(100*PeV);
      G4AntiOmegaMinusInelasticProcess* theAntiOmegaMinusInelasticProcess = new
G4AntiOmegaMinusInelasticProcess;
      theAntiOmegaMinusInelasticProcess->RegisterMe(theLEAntiOmegaMinusInelasticModel);
      theAntiOmegaMinusInelasticProcess->RegisterMe(theHEAntiOmegaMinusInelasticModel);
      //      theAntiOmegaMinusInelasticProcess->RegisterMe(theHEModel); //it
was missing before
      pmanager->AddDiscreteProcess(theAntiOmegaMinusInelasticProcess);

    }
  }
}

 old without using physics lists
-------------------------------------------------------- */

#endif

///////////////////////////////////////////////////////////////////////////
#include "G4OpAbsorption.hh"
#ifdef G4BOUNDARY_COMPILE
#include "G4OpBoundaryProcess.hh"
#endif
#ifdef G4ENABLE_MIE
#include "G4OpMie.hh"
#endif

void KM3Physics::ConstructOP() {
  theCerenkovProcess = new KM3Cherenkov("KM3Cherenkov");
  G4OpAbsorption *theAbsorptionProcess = new G4OpAbsorption();

  theCerenkovProcess->DumpPhysicsTable();
  // theScintillationProcess->DumpPhysicsTable();
  // theAbsorptionProcess->DumpPhysicsTable();

  theCerenkovProcess->SetVerboseLevel(0);
  theCerenkovProcess->SetDetector(aDetector);
  theAbsorptionProcess->SetVerboseLevel(0);

  G4int MaxNumPhotons = -30000;

  theCerenkovProcess->SetTrackSecondariesFirst(true);
  theCerenkovProcess->SetMaxNumPhotonsPerStep(MaxNumPhotons);

#ifdef G4BOUNDARY_COMPILE
  G4OpBoundaryProcess *theBoundaryProcess = new G4OpBoundaryProcess();
  theBoundaryProcess->SetVerboseLevel(0);
  G4OpticalSurfaceModel themodel = unified;
  theBoundaryProcess->SetModel(themodel);
#endif

#ifdef G4ENABLE_MIE
  G4OpMie *theMieProcess = new G4OpMie();
  theMieProcess->SetVerboseLevel(0);
#endif

  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition *particle = theParticleIterator->value();
    G4ProcessManager *pmanager = particle->GetProcessManager();
    if (pmanager) {
      G4String particleName = particle->GetParticleName();
      if (theCerenkovProcess->IsApplicable(*particle)) {
        pmanager->AddProcess(theCerenkovProcess);
        pmanager->SetProcessOrdering(theCerenkovProcess, idxPostStep);
      }
      if (particleName == "opticalphoton") {
        pmanager->AddDiscreteProcess(theAbsorptionProcess);
#ifdef G4BOUNDARY_COMPILE
        pmanager->AddDiscreteProcess(theBoundaryProcess);
#endif
#ifdef G4ENABLE_MIE
        pmanager->AddDiscreteProcess(theMieProcess);
#endif
      }
    }
  }
  // set ordering for AtRestDoIt
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//#include "G4Decay.hh"
// Decays ///////////////////////////////////////////////////////////////////
#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"

void KM3Physics::ConstructGeneral() {

  // Add Decay Process
  G4Decay *theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition *particle = theParticleIterator->value();
    G4ProcessManager *pmanager = particle->GetProcessManager();
    if (pmanager) {
      if (theDecayProcess->IsApplicable(*particle)) {
        // set ordering for PostStepDoIt and AtRestDoIt. Muons do not decay but
        // are captured (dense medium). Only for Hadronic runs
        if (particle->GetParticleName() == "mu-") {
#ifndef G4HADRONIC_COMPILE
          pmanager->AddProcess(theDecayProcess);
          pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
          pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
#endif
          ;
        } else {
          pmanager->AddProcess(theDecayProcess);
          pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
#ifdef G4HADRONIC_COMPILE
          if ((particle->GetParticleName() != "pi-") &&
              (particle->GetParticleName() != "kaon-"))
#endif
            pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
        }
      }
    }
  }

  // Declare radioactive decay to the GenericIon in the IonTable.
  const G4IonTable *theIonTable =
      G4ParticleTable::GetParticleTable()->GetIonTable();
  G4RadioactiveDecay *theRadioactiveDecay = new G4RadioactiveDecay();

  for (G4int i = 0; i < theIonTable->Entries(); i++) {
    G4String particleName = theIonTable->GetParticle(i)->GetParticleName();
    G4String particleType = theIonTable->GetParticle(i)->GetParticleType();
    // G4cout << "*********************************"<< particleName << G4endl;
    if (particleName == "GenericIon") {
      G4ProcessManager *pmanager =
          theIonTable->GetParticle(i)->GetProcessManager();
      pmanager->SetVerboseLevel(0);
      pmanager->AddProcess(theRadioactiveDecay);
      pmanager->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
      pmanager->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
    }
  }
}
// THIS IS END OF TEST

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//--start of apostolis parametrization----
void KM3Physics::AddParameterisation() {
  G4FastSimulationManagerProcess *theFastSimulationManagerProcess =
      new G4FastSimulationManagerProcess();

  //  theParticleIterator->reset();
  //  while( (*theParticleIterator)() ){
  //    G4ParticleDefinition* particle = theParticleIterator->value();
  //    if(!particle->IsShortLived()){         //ShortLived Particles have not
  //    process manager (they are not transported but decay immediately)
  //      G4ProcessManager* pmanager = particle->GetProcessManager();
  // both postStep and alongStep action are required: because
  // of the use of ghost volumes. If no ghost, the postStep
  // is sufficient.
  //      pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);
  //    }
  //  }

  // change and assign fast parametrization only to electron,positron and gamma
  G4ProcessManager *pmanager = G4Electron::Electron()->GetProcessManager();
  pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);

  pmanager = G4Positron::Positron()->GetProcessManager();
  pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);

  pmanager = G4Gamma::Gamma()->GetProcessManager();
  pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);

#ifdef G4HADRONIC_COMPILE
  // also add fast simulation to some hadrons

  pmanager = G4PionPlus::PionPlus()->GetProcessManager();
  pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);

  pmanager = G4PionMinus::PionMinus()->GetProcessManager();
  pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);

  pmanager = G4KaonPlus::KaonPlus()->GetProcessManager();
  pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);

  pmanager = G4KaonMinus::KaonMinus()->GetProcessManager();
  pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);

  pmanager = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();
  pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);

  pmanager = G4Proton::Proton()->GetProcessManager();
  pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);

  pmanager = G4AntiProton::AntiProton()->GetProcessManager();
  pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);

  pmanager = G4Neutron::Neutron()->GetProcessManager();
  pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);

  pmanager = G4AntiNeutron::AntiNeutron()->GetProcessManager();
  pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);

#endif
}
//--end of apostolis parametrization--------------

void KM3Physics::SetCuts() {
  // G4VUserPhysicsList::SetCutsWithDefault method sets
  // the default cut value for all particle types
  //
  SetCutsWithDefault();
  /*
  SetCutValue(defaultCutValue,"gamma");
  SetCutValue(0.6*mm,"e-");
  SetCutValue(0.6*mm,"e+");
  SetCutValue(defaultCutValue,"mu-");
  SetCutValue(defaultCutValue,"mu+");
  SetCutValue(defaultCutValue,"proton");
  SetCutValue(defaultCutValue,"anti_proton");
  SetCutValue(100*m,"opticalphoton");
  SetCutValueForOthers(defaultCutValue);
  */
  /*
  SetCutValue(1.1*m,"gamma");
  SetCutValue(0.6*mm,"e-");
  SetCutValue(0.6*mm,"e+");
  SetCutValue(7.0*m,"mu-");
  SetCutValue(7.0*m,"mu+");
  SetCutValue(defaultCutValue,"proton");
  SetCutValue(defaultCutValue,"anti_proton");
  SetCutValueForOthers(defaultCutValue);
  */

  printf("dfdfdfdfdfdfdfdf\n");
  DumpCutValuesTable(5); // b quark
  DumpCutValuesTable(22); // gamma
  DumpCutValuesTable(-11); // e+
  DumpCutValuesTable(11); // e-
  DumpCutValuesTable(-13); // mu+
  DumpCutValuesTable(13); // mu-
  DumpCutValuesTable(0); // opticalphoton
  printf("sssssssssssssssss\n");

  //  DumpCutValues();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
