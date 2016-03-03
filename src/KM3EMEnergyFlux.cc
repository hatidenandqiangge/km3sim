#include "KM3EMEnergyFlux.hh"
#include "Randomize.hh"
#include "CLHEP/Random/RandGamma.h"
#include "CLHEP/Random/RandPoisson.h"

KM3EMEnergyFlux::KM3EMEnergyFlux(char *infileParam, G4double QEmax,
                                 G4double TotCathodArea, G4int NEner,
                                 G4double MBstRat) {
  NEnergies = NEner;
  MaxBoostRatio = MBstRat;
  std::ifstream infile(infileParam, std::ios::in | std::ios::binary);
  keepEnergies = new std::vector<KM3EMDistanceFlux *>;
  keepEnergies->reserve(NEnergies);
  EnergyMin = 1.0e20;
  EnergyMax = -1.0e20;
  for (G4int i = 0; i < NEnergies; i++) {
    KM3EMDistanceFlux *aDistanceFlux = new KM3EMDistanceFlux(infile);
    if (EnergyMin > aDistanceFlux->GiveEnergy())
      EnergyMin = aDistanceFlux->GiveEnergy();
    if (EnergyMax < aDistanceFlux->GiveEnergy())
      EnergyMax = aDistanceFlux->GiveEnergy();
    keepEnergies->push_back(aDistanceFlux);
  }
  infile.close();
  RatioThis = QEmax * TotCathodArea;
}
KM3EMEnergyFlux::~KM3EMEnergyFlux() {
  if (keepEnergies != NULL) {
    for (G4int i = 0; i < NEnergies; i++) delete (*keepEnergies)[i];
    keepEnergies->clear();
    delete keepEnergies;
    keepEnergies = NULL;
  }
}

G4bool KM3EMEnergyFlux::ModelTrigger(G4double TheE) {
  if (log10(TheE) < EnergyMin || log10(TheE / MaxBoostRatio) > EnergyMax)
    return false;  // linear extrapolation is used for up to MaxBoostRatio times
                   // the maximum energy
  else
    return true;
}

// we have commented out the use of RMS sinc it is calculated wrongly in the
// parametrization section
// (see KM3SD and KM3EventAction)
// so we use strictly poisson distribution until it is fixed
void KM3EMEnergyFlux::FindBins(G4double energyin, G4double distancein,
                               G4double anglein) {
  G4double TheE = log10(energyin);
  G4int i;
  G4double BoostRatio;
  if (TheE < EnergyMax) {
    for (i = 1; i < NEnergies; i++) {
      if (TheE < (*keepEnergies)[i]->GiveEnergy()) break;
    }
    ibin2 = i;
    BoostRatio = 1.0;
  } else {
    ibin2 = NEnergies - 1;
    BoostRatio = energyin / pow(10.0, EnergyMax);
    TheE = EnergyMax;
  }
  ibin1 = ibin2 - 1;

  (*keepEnergies)[ibin1]->FindBins(distancein, anglein);
  (*keepEnergies)[ibin2]->FindBins(distancein, anglein);

  ratio = (TheE - (*keepEnergies)[ibin1]->GiveEnergy()) /
          ((*keepEnergies)[ibin2]->GiveEnergy() -
           (*keepEnergies)[ibin1]->GiveEnergy());
  G4double Flux1 = (*keepEnergies)[ibin1]->GiveFlux();
  G4double Flux2 = (*keepEnergies)[ibin2]->GiveFlux();
  Flux = Flux1 + ratio * (Flux2 - Flux1);

  G4double dFlux02 = 1.0 - exp(Flux - Flux2);
  G4double dFlux12 = 1.0 - exp(Flux1 - Flux2);
  ratio = dFlux02 / dFlux12;
  if (ratio < 0)
    ratio = 0.0;
  else if (ratio > 1.0)
    ratio = 1.0;

  Flux = BoostRatio * RatioThis * exp(Flux);

  // tempoR  FluxRMS=(*keepEnergies)[ibin1]->GiveFluxRMS() +
  // ratio*((*keepEnergies)[ibin2]->GiveFluxRMS() -
  // (*keepEnergies)[ibin1]->GiveFluxRMS() );
  // tempoR  FluxRMS=BoostRatio*RatioThis*exp(FluxRMS);
  // next sample according to the Polya (also called Pascal or negative
  // binomial) distribution
  // find the number of samples according to Flux and FluxRMS
  // tempoR  G4double deviation=FluxRMS-Flux;// It is variance-mean.this should
  // be >0 for polya, =0 for poisson.
  // tempoR  if(deviation < 0.001){ //this case is for poisson (the limiting
  // case of polya as n->oo is the poisson).
  NumberOfSamples = (G4int)CLHEP::RandPoisson::shoot(Flux);
  // tempoR  }
  // tempoR  else{
  // tempoR    G4double p=Flux/FluxRMS; //this is mean/variance
  // tempoR    G4double n=Flux*Flux/deviation;
  // tempoR    G4double X = CLHEP::RandGamma::shoot(n,1.0);
  // tempoR    NumberOfSamples = (G4int) CLHEP::RandPoisson::shoot(X*(1-p)/p);
  // tempoR  }
  //  printf("inside %le %le %d\n",meannum,rmsnum,numofphotons);
}

onePE KM3EMEnergyFlux::GetSamplePoint() {
  if (G4UniformRand() < ratio)
    return (*keepEnergies)[ibin1]->GetSamplePoint();
  else
    return (*keepEnergies)[ibin2]->GetSamplePoint();
}
