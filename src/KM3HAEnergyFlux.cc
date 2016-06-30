#include "KM3HAEnergyFlux.h"
#include "Randomize.hh"
#include "CLHEP/Random/RandGamma.h"
#include "CLHEP/Random/RandPoisson.h"

using CLHEP::TeV;
using CLHEP::GeV;
using CLHEP::meter;

KM3HAEnergyFlux::KM3HAEnergyFlux(char *infileParam, G4double QEmax,
                                 G4double TotCathodArea, G4double EneMin,
                                 G4double EneMax) {
  NPartsDists = 2;  // is the number of distributions kept (now it is 2, for 211
                    // (pi+) and 130 (KaonZeroLong))
  std::ifstream infile(infileParam, std::ios::in | std::ios::binary);
  keepEnergies = new std::vector<KM3EMDistanceFlux *>;
  keepEnergies->reserve(NPartsDists);
  EnergyMin = log10(EneMin);  // the energy range the distributions apply
  EnergyMax = log10(EneMax);
  for (G4int i = 0; i < NPartsDists; i++) {
    KM3EMDistanceFlux *aDistanceFlux = new KM3EMDistanceFlux(infile);
    keepEnergies->push_back(aDistanceFlux);
  }
  infile.close();
  RatioThis = QEmax * TotCathodArea;
  // next is photon output of HA particles with arbitrary scaling
  // they are for particles (-2212, -2112, 130, +-211, +-321, 2112, 2212) in
  // this order
  // it has been fitted (flux/(energy/TeV) vs log10(energy/TeV) with a poly 4th
  // order
  //-2212
  ParticleEnergyScaleFactors[0][0] = 0.12710;
  ParticleEnergyScaleFactors[0][1] = 0.14127e-01;
  ParticleEnergyScaleFactors[0][2] = -0.17273e-02;
  ParticleEnergyScaleFactors[0][3] = -0.67207e-03;
  ParticleEnergyScaleFactors[0][4] = 0.16668e-03;
  //-2112 not too good
  ParticleEnergyScaleFactors[1][0] = 0.12956;
  ParticleEnergyScaleFactors[1][1] = 0.15708e-01;
  ParticleEnergyScaleFactors[1][2] = -0.47459e-02;
  ParticleEnergyScaleFactors[1][3] = -0.93359e-03;
  ParticleEnergyScaleFactors[1][4] = 0.73662e-03;
  // 130
  ParticleEnergyScaleFactors[2][0] = 0.13744;
  ParticleEnergyScaleFactors[2][1] = 0.10331e-01;
  ParticleEnergyScaleFactors[2][2] = -0.35095e-02;
  ParticleEnergyScaleFactors[2][3] = 0.15392e-02;
  ParticleEnergyScaleFactors[2][4] = -0.26142e-03;
  //+-211
  ParticleEnergyScaleFactors[3][0] = 0.13656;
  ParticleEnergyScaleFactors[3][1] = 0.11316e-01;
  ParticleEnergyScaleFactors[3][2] = -0.14285e-02;
  ParticleEnergyScaleFactors[3][3] = 0.58399e-03;
  ParticleEnergyScaleFactors[3][4] = -0.46483e-03;
  //+-321
  ParticleEnergyScaleFactors[4][0] = 0.13230;
  ParticleEnergyScaleFactors[4][1] = 0.12236e-01;
  ParticleEnergyScaleFactors[4][2] = -0.14947e-02;
  ParticleEnergyScaleFactors[4][3] = 0.98115e-03;
  ParticleEnergyScaleFactors[4][4] = -0.55477e-03;
  // 2112
  ParticleEnergyScaleFactors[5][0] = 0.13080;
  ParticleEnergyScaleFactors[5][1] = 0.13472e-01;
  ParticleEnergyScaleFactors[5][2] = -0.21558e-02;
  ParticleEnergyScaleFactors[5][3] = 0.18663e-02;
  ParticleEnergyScaleFactors[5][4] = -0.82554e-03;
  // 2212
  ParticleEnergyScaleFactors[6][0] = 0.13130;
  ParticleEnergyScaleFactors[6][1] = 0.12880e-01;
  ParticleEnergyScaleFactors[6][2] = -0.28592e-02;
  ParticleEnergyScaleFactors[6][3] = 0.19248e-02;
  ParticleEnergyScaleFactors[6][4] = -0.61522e-03;
}

KM3HAEnergyFlux::~KM3HAEnergyFlux() {
  if (keepEnergies != NULL) {
    for (G4int i = 0; i < NPartsDists; i++) delete (*keepEnergies)[i];
    keepEnergies->clear();
    delete keepEnergies;
    keepEnergies = NULL;
  }
}

G4double KM3HAEnergyFlux::Rescale(G4int idbeam, G4double energyin,
                                  G4double energydis) {
  G4int irbin;
  // they are for particles (-2212, -2112, 130, +-211, +-321, 2112, 2212) in
  // this order
  if (idbeam == -2212)
    irbin = 0;
  else if (idbeam == -2112)
    irbin = 1;
  else if (idbeam == 130)
    irbin = 2;
  else if (idbeam == 211 || idbeam == -211)
    irbin = 3;
  else if (idbeam == 321 || idbeam == -321)
    irbin = 4;
  else if (idbeam == 2112)
    irbin = 5;
  else if (idbeam == 2212)
    irbin = 6;
  G4double eneinlog = log10(energyin / TeV);
  G4double enedislog = log10(energydis / TeV);
  G4double ydis =
      ParticleEnergyScaleFactors[irbin][0] +
      enedislog *
          (ParticleEnergyScaleFactors[irbin][1] +
           enedislog *
               (ParticleEnergyScaleFactors[irbin][2] +
                enedislog *
                    (ParticleEnergyScaleFactors[irbin][3] +
                     enedislog * (ParticleEnergyScaleFactors[irbin][4]))));
  G4double yin =
      ParticleEnergyScaleFactors[irbin][0] +
      eneinlog *
          (ParticleEnergyScaleFactors[irbin][1] +
           eneinlog *
               (ParticleEnergyScaleFactors[irbin][2] +
                eneinlog *
                    (ParticleEnergyScaleFactors[irbin][3] +
                     eneinlog * (ParticleEnergyScaleFactors[irbin][4]))));
  return (yin / ydis) * (energyin / energydis);
}

// we have commented out the use of RMS sinc it is calculated wrongly in the
// parametrization section
// (see KM3SD and KM3EventAction)
// so we use strictly poisson distribution until it is fixed
void KM3HAEnergyFlux::FindBins(G4int idbeam, G4double energyin,
                               G4double distancein, G4double anglein) {
  // we have 2 catogories of hadronic particles for the parametrization. Charged
  // ones and neutral ones
  // since they show very different angular profiles
  if (idbeam == 211 || idbeam == -211 || idbeam == 2212 || idbeam == -2212 ||
      idbeam == 321 || idbeam == -321)
    ibin = 0;
  else if (idbeam == 130 || idbeam == 2112 || idbeam == -2112)
    ibin = 1;

  (*keepEnergies)[ibin]->FindBins(distancein, anglein);

  G4double energydis = pow(10.0, (*keepEnergies)[ibin]->GiveEnergy());
  G4double BoostRatio = Rescale(idbeam, energyin, energydis);

  Flux = (*keepEnergies)[ibin]->GiveFlux();

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

onePE KM3HAEnergyFlux::GetSamplePoint() {
  return (*keepEnergies)[ibin]->GetSamplePoint();
}
