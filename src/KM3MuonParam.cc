#include "KM3MuonParam.hh"
#include <stdio.h>

KM3MuonParam::KM3MuonParam()
{
  FILE* infile;
  G4double energy,distance,prob0;
  G4int nbins;
  G4double tempoArray[10000];
  EnergyCutoff=250.0*GeV; //100GeV before
  if((infile=fopen("SimLengthProb","r")) ==NULL){
    G4Exception("Error open input Range-Energy Muons file\n","",FatalException,"");
  }
  else{
    MinLogEnergy=1.0e9;
    MaxLogEnergy=-1.0e9;
    while(fscanf(infile,"%lf %lf %lf %d\n",&energy,&distance,&prob0,&nbins) != EOF){
      energy*=GeV;
      for(G4int i=0 ; i< nbins ; i++)fscanf(infile,"%lf\n",&tempoArray[i]);
      //this code is for the redefinition of the lowest cutoff
      //first redefine the the prob0
      G4double numofentries=0;
      for(G4int i=0 ; i< nbins ; i++)numofentries += tempoArray[i];
      G4double numrejected=0;
      G4int iterat=1;
      G4double stopbin=energy*double(iterat+1)/double(nbins);
      while (stopbin <= EnergyCutoff){
	numrejected += tempoArray[iterat];
	tempoArray[iterat]=0.0;
	iterat++;
	stopbin=energy*double(iterat+1)/double(nbins);
      }
      G4double nlastrej=tempoArray[iterat]*(stopbin-EnergyCutoff)*double(nbins)/energy;
      numrejected += nlastrej;
      tempoArray[iterat] -= nlastrej;
      prob0=prob0+(1.0-prob0)*numrejected/numofentries;
      if(prob0<0.999){
	////////////////////////////////////////////////////
	PDFSList* aPDFSList=(PDFSList*)malloc(sizeof(PDFSList));
	aPDFSList->LogEnergy=log10(energy);
	aPDFSList->Distance=distance*m;
	aPDFSList->Prob0=prob0;
	aPDFSList->StopFirstBin=double(int((EnergyCutoff/energy)*double(nbins))+1)*energy/double(nbins);
	aPDFSList->thePDF = new CLHEP::RandGeneral(tempoArray,nbins);
	thePDFS.push_back(aPDFSList);
	if(aPDFSList->LogEnergy < MinLogEnergy)MinLogEnergy=aPDFSList->LogEnergy;
	if(aPDFSList->LogEnergy > MaxLogEnergy)MaxLogEnergy=aPDFSList->LogEnergy;
      }
    }
  }
  fclose(infile);
  MyTolerance=1.0e-3;
}

KM3MuonParam::~KM3MuonParam()
{
  for(size_t ien=0 ; ien<thePDFS.size() ; ien++){
    delete thePDFS[ien]->thePDF;
    free(thePDFS[ien]);
  }
  thePDFS.clear();
}

void KM3MuonParam::AddMuon(G4double energy, G4double distance)
{
  G4double logenergy=log10(energy);
  forMuon* aforMuon=(forMuon*)malloc(sizeof(forMuon));
  aforMuon->LogEnergy=logenergy;
  aforMuon->Distance=0.0;
  aforMuon->Prob0=1.0;
  aforMuon->thePDF=NULL;
  aforMuon->iscapable=false;
  if(logenergy > MaxLogEnergy){
    aforMuon->iscapable=true;
  }
  if(logenergy>MinLogEnergy && logenergy<MaxLogEnergy){
    aforMuon->iscapable=true;
    //first find the closest energy in logarithm
    G4double mindist=1.e9;
    G4double enethis;
    for(size_t ip=0 ; ip<thePDFS.size() ; ip++){
      if(fabs(logenergy-thePDFS[ip]->LogEnergy) < mindist){
	mindist=fabs(logenergy-thePDFS[ip]->LogEnergy);
	enethis=thePDFS[ip]->LogEnergy;
      }
    }
    //next for the closest energy find the closest spatial distance
    mindist=1.e9;
    size_t ipthis;
    for(size_t ip=0 ; ip<thePDFS.size() ; ip++){
      if(fabs(enethis-thePDFS[ip]->LogEnergy) < MyTolerance){
	if(fabs(distance-thePDFS[ip]->Distance) < mindist){
	  mindist=fabs(distance-thePDFS[ip]->Distance);
	  ipthis=ip; 
	}
      }
    }
    if(distance-thePDFS[ipthis]->Distance > 100.0*m)aforMuon->iscapable=false;
    else if(distance-thePDFS[ipthis]->Distance > -100.0*m){
      aforMuon->thePDF=thePDFS[ipthis]->thePDF;
      aforMuon->Distance=thePDFS[ipthis]->Distance;
      aforMuon->Prob0=thePDFS[ipthis]->Prob0;
      aforMuon->StopFirstBin=thePDFS[ipthis]->StopFirstBin;
    }
  }
  theDistributions.push_back(aforMuon);
}
void KM3MuonParam::Initialize(void)
{
  isEventOK=true;
  for(size_t ip=0 ; ip<theDistributions.size() ; ip++){
    if(theDistributions[ip]->iscapable && theDistributions[ip]->thePDF==NULL){
      isEventOK=false;
      break;
    }
  }
  if(!isEventOK){
    for(size_t ip=0 ; ip<theDistributions.size() ; ip++){
      theDistributions[ip]->Prob0=0.0;
      theDistributions[ip]->Distance=0.0;
    }
  }
}

void KM3MuonParam::Finalize(void)
{
  for(size_t ip=0 ; ip<theDistributions.size() ; ip++)
    free(theDistributions[ip]);
  theDistributions.clear();
}

G4double KM3MuonParam::GetDistance(G4int idmuon)
{return theDistributions[idmuon]->Distance;}

G4double KM3MuonParam::GetEnergy(G4int idmuon)
{
  if(!isEventOK)
    return pow(10.0,theDistributions[idmuon]->LogEnergy);
  else{
    if(G4UniformRand()<theDistributions[idmuon]->Prob0)return 0.0;
    else {
      G4double thefrac=theDistributions[idmuon]->thePDF->shoot();
      G4double theenergy=thefrac*pow(10.0,theDistributions[idmuon]->LogEnergy);
      if(theenergy<EnergyCutoff)
	theenergy=EnergyCutoff+(theDistributions[idmuon]->StopFirstBin-EnergyCutoff)*G4UniformRand();
      return theenergy;
    }
  }
}

G4bool KM3MuonParam::IsCapable(G4int idmuon)
{return theDistributions[idmuon]->iscapable;}

G4double KM3MuonParam::GetWeight(void)
{
  G4double allprob0=1.0;
  for(size_t ip=0 ; ip<theDistributions.size() ; ip++)
    allprob0 *= theDistributions[ip]->Prob0;
  return 1.0-allprob0;
}
