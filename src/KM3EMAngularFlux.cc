#include "KM3EMAngularFlux.hh"
#include "Randomize.hh"


KM3EMAngularFlux::KM3EMAngularFlux(std::ifstream& infile,bool& ok,bool FineBin)
{
  VertexSolidAngleBins=51;
  if(FineBin){
    VertexSolidAngleBins=71;
  }
  keepAngles = new std::vector<KM3EMTimePointDis*>;
  keepAngles->reserve(VertexSolidAngleBins);
  char valC[4];
  infile.read(valC,4);
  Distance=double(*(float*)valC);
  Distance *= meter;
  G4int count=0;
  G4int countall=0;
  for(G4int i=0 ; i<VertexSolidAngleBins ; i++){
    bool oka;
    KM3EMTimePointDis* aTimePointDis= new KM3EMTimePointDis(infile,oka);
    keepAngles->push_back(aTimePointDis);
    if(oka)count=0;
    else count++;
    if(count > 1)countall++;
    if(!oka) G4cout <<"Null for Distance " <<Distance/meter<<" and angle "<< aTimePointDis->GiveAngle() <<G4endl;
  }
  if(countall>4)ok=false;
  else ok=true;
  if(ok)IsThisValid=true;
  else{
    IsThisValid=false;
    for(G4int i=0 ; i<VertexSolidAngleBins ; i++)delete (*keepAngles)[i];
    keepAngles->clear();
    delete keepAngles;
    keepAngles=NULL;
  }
}
KM3EMAngularFlux::~KM3EMAngularFlux()
{
  if(keepAngles!=NULL){
    for(G4int i=0 ; i<VertexSolidAngleBins ; i++)delete (*keepAngles)[i];
    keepAngles->clear();
    delete keepAngles;
    keepAngles=NULL;
  }
}
void KM3EMAngularFlux::FindBins(G4double anglein)
{
  if(!IsThisValid)G4Exception("Error sampling angle for null distribution\n","",FatalException,"");
  G4int i;
  for(i=0 ; i<VertexSolidAngleBins ; i++){
    if((*keepAngles)[i]->IsValid()){
      if(anglein < (*keepAngles)[i]->GiveAngle())break;
    }
  }
  ibin2=i;
  if(ibin2==VertexSolidAngleBins){
    for(i=VertexSolidAngleBins-1 ; i>=0 ; i--){
      if((*keepAngles)[i]->IsValid())break;
    }
    ibin2=i;
  }
  for(i=ibin2-1 ; i>=0 ; i--){
    if((*keepAngles)[i]->IsValid())break;
  }
  ibin1=i;
  if(ibin1 == -1){
    ibin1=ibin2;
    for(i=ibin1+1 ; i<VertexSolidAngleBins ; i++){
      if((*keepAngles)[i]->IsValid())break;
    }
    ibin2=i;
  }
  if(ibin2 ==0){
    ibin1=ibin2;
    for(i=ibin1+1 ; i<VertexSolidAngleBins ; i++){
      if((*keepAngles)[i]->IsValid())break;
    }
    ibin2=i;
  }
  G4double Angle1=(*keepAngles)[ibin1]->GiveAngle();
  G4double Angle2=(*keepAngles)[ibin2]->GiveAngle();
  ratio=(anglein-Angle1)/(Angle2 - Angle1);
  G4double Flux1=(*keepAngles)[ibin1]->GiveFlux();
  G4double Flux2=(*keepAngles)[ibin2]->GiveFlux();
  Flux=Flux1 + ratio*(Flux2 - Flux1);
//tempoR  FluxRMS=(*keepAngles)[ibin1]->GiveFluxRMS() + ratio*((*keepAngles)[ibin2]->GiveFluxRMS() - (*keepAngles)[ibin1]->GiveFluxRMS() );
  G4double dFlux02=1.0-exp(Flux -Flux2);
  G4double dFlux12=1.0-exp(Flux1-Flux2);
  ratio = dFlux02/dFlux12;
  if(ratio<0)ratio=0.0;
  else if(ratio>1.0)ratio=1.0;
}

onePE KM3EMAngularFlux::GetSamplePoint()
{
  if(G4UniformRand()<ratio)return (*keepAngles)[ibin1]->GetSamplePoint();
  else return (*keepAngles)[ibin2]->GetSamplePoint();
}

