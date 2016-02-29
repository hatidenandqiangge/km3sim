#include "KM3EMDirectFlux.hh"
#include "Randomize.hh"


KM3EMDirectFlux::KM3EMDirectFlux(char* infileParam,G4double TotCathodArea)
{
  VertexDistanceBins=40;
  std::ifstream infile(infileParam, std::ios::in | std::ios::binary);
  //keep2013  infile.seekg(std::streampos(439203908));
  infile.seekg(std::streampos(540323872)); //it is 67540484*8 where the first number is size/energy and the second the number of energies for e-/e+/gamma parametrization
  keepDistances = new std::vector<KM3EMAngularFlux*>;
  keepDistances->reserve(VertexDistanceBins);
  char valC[4];
  infile.read(valC,4);
  G4double Energy=double(*(float*)valC);
  for(G4int i=0 ; i<VertexDistanceBins ; i++){
    bool oka;
    KM3EMAngularFlux* aAngularFlux= new KM3EMAngularFlux(infile,oka,true); //true is for fine binning
    keepDistances->push_back(aAngularFlux);
    if(!oka) G4cout <<"Null for Energy " <<Energy<<" and distance "<< aAngularFlux->GiveDistance()/meter <<G4endl;
  }
  infile.close();
  RatioThis=TotCathodArea;
}
KM3EMDirectFlux::~KM3EMDirectFlux()
{
  if(keepDistances!=NULL){
    for(G4int i=0 ; i<VertexDistanceBins ; i++)delete (*keepDistances)[i];
    keepDistances->clear();
    delete keepDistances;
    keepDistances=NULL;
  }
}
void KM3EMDirectFlux::FindBins(G4double MeanNumPhotons,G4double distancein, G4double anglein)
{
  G4int i;
  for(i=0 ; i<VertexDistanceBins ; i++){
    if((*keepDistances)[i]->IsValid()){
      if(distancein < (*keepDistances)[i]->GiveDistance())break;
    }
  }
  ibin2=i;
  if(ibin2==VertexDistanceBins){
    for(i=VertexDistanceBins-1 ; i>=0 ; i--){
      if((*keepDistances)[i]->IsValid())break;
    }
    ibin2=i;
  }
  for(i=ibin2-1 ; i>=0 ; i--){
    if((*keepDistances)[i]->IsValid())break;
  }
  ibin1=i;
  if(ibin1 == -1){
    ibin1=ibin2;
    for(i=ibin1+1 ; i<VertexDistanceBins ; i++){
      if((*keepDistances)[i]->IsValid())break;
    }
    ibin2=i;
  }
  if(ibin2 ==0){
    ibin1=ibin2;
    for(i=ibin1+1 ; i<VertexDistanceBins ; i++){
      if((*keepDistances)[i]->IsValid())break;
    }
    ibin2=i;
  }
  (*keepDistances)[ibin1]->FindBins(anglein);
  (*keepDistances)[ibin2]->FindBins(anglein);

  //new d^2*F interpolation
  G4double Distance1=(*keepDistances)[ibin1]->GiveDistance();
  G4double Distance2=(*keepDistances)[ibin2]->GiveDistance();
  ratio=(distancein-Distance1)/(Distance2 - Distance1);
  G4double Flux1=(*keepDistances)[ibin1]->GiveFlux();
  G4double Flux2=(*keepDistances)[ibin2]->GiveFlux();
  Flux=Flux1+ratio*(Flux2-Flux1);
  Flux+=2*(log(Distance1/distancein)+ratio*log(Distance2/Distance1));
  //new d^2F interpolation
  G4double dFlux02=1.0-exp(Flux -Flux2);
  G4double dFlux12=1.0-exp(Flux1-Flux2);
  ratio = dFlux02/dFlux12;
  if(ratio<0)ratio=0.0;
  else if(ratio>1.0)ratio=1.0;
  Flux=RatioThis*MeanNumPhotons*exp(Flux);
  NumberOfSamples = (G4int) CLHEP::RandPoisson::shoot(Flux);
}

onePE KM3EMDirectFlux::GetSamplePoint()
{
  if(G4UniformRand()<ratio)return(*keepDistances)[ibin1]->GetSamplePoint();
  else return (*keepDistances)[ibin2]->GetSamplePoint();
}

