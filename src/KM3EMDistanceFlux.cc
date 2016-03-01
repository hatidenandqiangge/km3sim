#include "KM3EMDistanceFlux.hh"
#include "Randomize.hh"

KM3EMDistanceFlux::KM3EMDistanceFlux(std::ifstream &infile) {
  VertexDistanceBins = 40;
  keepDistances = new std::vector<KM3EMAngularFlux *>;
  keepDistances->reserve(VertexDistanceBins);
  char valC[4];
  infile.read(valC, 4);
  Energy = log10(GeV * double(*(float *)valC));
  for (G4int i = 0; i < VertexDistanceBins; i++) {
    bool oka;
    KM3EMAngularFlux *aAngularFlux =
        new KM3EMAngularFlux(infile, oka, false); // false is for fine binning
    keepDistances->push_back(aAngularFlux);
    if (!oka)
      G4cout << "Null for Energy " << Energy << " and distance "
             << aAngularFlux->GiveDistance() / meter << G4endl;
  }
}
KM3EMDistanceFlux::~KM3EMDistanceFlux() {
  if (keepDistances != NULL) {
    for (G4int i = 0; i < VertexDistanceBins; i++)
      delete (*keepDistances)[i];
    keepDistances->clear();
    delete keepDistances;
    keepDistances = NULL;
  }
}
void KM3EMDistanceFlux::FindBins(G4double distancein, G4double anglein) {
  G4int i;
  for (i = 0; i < VertexDistanceBins; i++) {
    if ((*keepDistances)[i]->IsValid()) {
      if (distancein < (*keepDistances)[i]->GiveDistance())
        break;
    }
  }
  ibin2 = i;
  if (ibin2 == VertexDistanceBins) {
    for (i = VertexDistanceBins - 1; i >= 0; i--) {
      if ((*keepDistances)[i]->IsValid())
        break;
    }
    ibin2 = i;
  }
  for (i = ibin2 - 1; i >= 0; i--) {
    if ((*keepDistances)[i]->IsValid())
      break;
  }
  ibin1 = i;
  if (ibin1 == -1) {
    ibin1 = ibin2;
    for (i = ibin1 + 1; i < VertexDistanceBins; i++) {
      if ((*keepDistances)[i]->IsValid())
        break;
    }
    ibin2 = i;
  }
  if (ibin2 == 0) {
    ibin1 = ibin2;
    for (i = ibin1 + 1; i < VertexDistanceBins; i++) {
      if ((*keepDistances)[i]->IsValid())
        break;
    }
    ibin2 = i;
  }
  (*keepDistances)[ibin1]->FindBins(anglein);
  (*keepDistances)[ibin2]->FindBins(anglein);

  // new d^2*F interpolation
  G4double Distance1 = (*keepDistances)[ibin1]->GiveDistance();
  G4double Distance2 = (*keepDistances)[ibin2]->GiveDistance();
  ratio = (distancein - Distance1) / (Distance2 - Distance1);
  G4double Flux1 = (*keepDistances)[ibin1]->GiveFlux();
  G4double Flux2 = (*keepDistances)[ibin2]->GiveFlux();
  Flux = Flux1 + ratio * (Flux2 - Flux1);
  Flux +=
      2 * (log(Distance1 / distancein) + ratio * log(Distance2 / Distance1));
  // new d^2F interpolation
  // tempoR  FluxRMS=(*keepDistances)[ibin1]->GiveFluxRMS() +
  // ratio*((*keepDistances)[ibin2]->GiveFluxRMS() -
  // (*keepDistances)[ibin1]->GiveFluxRMS() );
  G4double dFlux02 = 1.0 - exp(Flux - Flux2);
  G4double dFlux12 = 1.0 - exp(Flux1 - Flux2);
  ratio = dFlux02 / dFlux12;
  if (ratio < 0)
    ratio = 0.0;
  else if (ratio > 1.0)
    ratio = 1.0;
}

onePE KM3EMDistanceFlux::GetSamplePoint() {
  if (G4UniformRand() < ratio)
    return (*keepDistances)[ibin1]->GetSamplePoint();
  else
    return (*keepDistances)[ibin2]->GetSamplePoint();
}
