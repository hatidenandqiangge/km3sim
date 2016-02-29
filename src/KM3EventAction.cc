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
// $Id: G4UserEventAction.cc,v 1.7 2006/06/29 18:10:31 gunter Exp $
// GEANT4 tag $Name: geant4-08-01-patch-01 $
//

#include "KM3EventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"


void KM3EventAction::BeginOfEventAction(const G4Event*)
{
  if(!(G4ParticleTable::GetParticleTable()->GetReadiness()))
 {
   G4String msg;
   msg =  " You are instantiating G4UserEventAction BEFORE your\n";
   msg += "G4VUserPhysicsList is instantiated and assigned to G4RunManager.\n";
   msg += " Such an instantiation is prohibited by Geant4 version 8.0. To fix this problem,\n";
   msg += "please make sure that your main() instantiates G4VUserPhysicsList AND\n";
   msg += "set it to G4RunManager before instantiating other user action classes\n";
   msg += "such as G4UserEventAction.";
   G4Exception("G4UserEventAction::G4UserEventAction()",
              "Event0001",FatalException,msg);
 }
  centerPre.clear();
  centerPost.clear();
  enterPre.clear();
  enterPost.clear();
  leavePre.clear();
  leavePost.clear();
  centerMomentum.clear();
  enterMomentum.clear();
  leaveMomentum.clear();
  centerPosition.clear();
  enterPosition.clear();
  leavePosition.clear();
  centerTime.clear();
  enterTime.clear();
  leaveTime.clear();
  stopPosition.clear();
  stopTime.clear();
#ifdef G4MYMUON_KEEPENERGY
  EnergyAtPosition.clear();
#endif

  G4ThreeVector vzero(0.0,0.0,0.0);
  G4double izero;
  izero=0.0;
  for(G4int ip=0 ; ip<numofMuons ; ip++){
    centerPre.push_back(vzero);
    centerPost.push_back(vzero);
    enterPre.push_back(vzero);
    enterPost.push_back(vzero);
    leavePre.push_back(vzero);
    leavePost.push_back(vzero);
    centerMomentum.push_back(izero);
    enterMomentum.push_back(izero);
    leaveMomentum.push_back(izero);
    centerPosition.push_back(vzero);
    enterPosition.push_back(vzero);
    leavePosition.push_back(vzero);
    centerTime.push_back(izero);
    enterTime.push_back(izero);
    leaveTime.push_back(izero);
    stopPosition.push_back(vzero);
    stopTime.push_back(izero);
  }
  if(useANTARESformat)TheEVTtoWrite->ReadEvent();
#if defined(G4MYEM_PARAMETERIZATION) || defined(G4MYHA_PARAMETERIZATION) //newha
  G4int TotalNumberOfCathods=MyStDetector->allCathods->GetNumberOfCathods();
  G4bool FineBin=false;
  G4int VertexSolidAngleBins=51;
  if(MyGenerator->ParamEnergy == 0.0){
    FineBin=true;
    VertexSolidAngleBins=71;
  }
  for(G4int icath=0 ; icath<TotalNumberOfCathods ; icath++){
    for(G4int cang1bin=0 ; cang1bin<VertexSolidAngleBins ; cang1bin++){
      G4int distbin=icath;
      G4int ibin_d1=cang1bin+VertexSolidAngleBins*distbin;
      (*myPhotonsNumber)[ibin_d1]=0.0;
    }
  }
#endif
}


void KM3EventAction::EndOfEventAction(const G4Event*)
{
  //write the momentums, positions and times to out file
  G4ThreeVector Momentum;
  G4ThreeVector vzero(0.0,0.0,0.0);
  G4double izero=0.0;
#ifndef G4MYK40_PARAMETERIZATION
  if(!useANTARESformat)fprintf(outfile,"%d\n",numofMuons);
#endif
  for(G4int ip=0 ; ip<numofMuons ; ip++){
    //record entering position//////
    if( (enterPre[ip] != vzero) && (enterPost[ip] != vzero)){
      Momentum=((enterMomentum[ip]))*(enterPost[ip]-enterPre[ip])/(enterPost[ip]-enterPre[ip]).mag();
      if(!useANTARESformat)fprintf(outfile,"%.6e %.6e %.6e %.6e %.6e %.6e %.7e\n",
	      (enterPosition[ip])[0]/cm,(enterPosition[ip])[1]/cm,(enterPosition[ip])[2]/cm,
	      Momentum[0]/GeV,Momentum[1]/GeV,Momentum[2]/GeV,
	      (enterTime[ip])/second);
      else {
	Momentum=Momentum.unit();
	TheEVTtoWrite->AddMuonPositionInfo(MuonIds[ip],-1,
					   (enterPosition[ip])[0]/m,(enterPosition[ip])[1]/m,(enterPosition[ip])[2]/m,
					   Momentum[0],Momentum[1],Momentum[2],enterMomentum[ip]/GeV,
					   (enterTime[ip])/ns);
      }
    }
    else{
      if(!useANTARESformat)fprintf(outfile,"%.6e %.6e %.6e %.6e %.6e %.6e %.7e\n",
	      vzero[0],vzero[1],vzero[2],
	      vzero[0],vzero[1],vzero[2],
	      izero);
      else TheEVTtoWrite->AddMuonPositionInfo(MuonIds[ip],-1,0.,0.,0.,0.,0.,0.,0.,0.);
    }
    //record center position///////////
    if( (centerPre[ip] != vzero) && (centerPost[ip] != vzero)){
      Momentum=((centerMomentum[ip]))*(centerPost[ip]-centerPre[ip])/(centerPost[ip]-centerPre[ip]).mag();
      if(!useANTARESformat)fprintf(outfile,"%.6e %.6e %.6e %.6e %.6e %.6e %.7e\n",
	      (centerPosition[ip])[0]/cm,(centerPosition[ip])[1]/cm,(centerPosition[ip])[2]/cm,
	      Momentum[0]/GeV,Momentum[1]/GeV,Momentum[2]/GeV,
	      (centerTime[ip])/second);
      else {
	Momentum=Momentum.unit();
	TheEVTtoWrite->AddMuonPositionInfo(MuonIds[ip],0,
					   (centerPosition[ip])[0]/m,(centerPosition[ip])[1]/m,(centerPosition[ip])[2]/m,
					   Momentum[0],Momentum[1],Momentum[2],centerMomentum[ip]/GeV,
					   (centerTime[ip])/ns);
      }
    }
    else{
      if(!useANTARESformat)fprintf(outfile,"%.6e %.6e %.6e %.6e %.6e %.6e %.7e\n",
	      vzero[0],vzero[1],vzero[2],
	      vzero[0],vzero[1],vzero[2],
	      izero);
      else TheEVTtoWrite->AddMuonPositionInfo(MuonIds[ip],0,0.,0.,0.,0.,0.,0.,0.,0.);
    }
    //record leaving position//////////////////////
    if( (leavePre[ip] != vzero) && (leavePost[ip] != vzero)){
      Momentum=((leaveMomentum[ip]))*(leavePost[ip]-leavePre[ip])/(leavePost[ip]-leavePre[ip]).mag();
      if(!useANTARESformat)fprintf(outfile,"%.6e %.6e %.6e %.6e %.6e %.6e %.7e\n",
	      (leavePosition[ip])[0]/cm,(leavePosition[ip])[1]/cm,(leavePosition[ip])[2]/cm,
	      Momentum[0]/GeV,Momentum[1]/GeV,Momentum[2]/GeV,
	      (leaveTime[ip])/second);
      else {
	Momentum=Momentum.unit();
	TheEVTtoWrite->AddMuonPositionInfo(MuonIds[ip],1,
					   (leavePosition[ip])[0]/m,(leavePosition[ip])[1]/m,(leavePosition[ip])[2]/m,
					   Momentum[0],Momentum[1],Momentum[2],leaveMomentum[ip]/GeV,
					   (leaveTime[ip])/ns);
      }
    }
    else{
      if(!useANTARESformat)fprintf(outfile,"%.6e %.6e %.6e %.6e %.6e %.6e %.7e\n",
	      vzero[0],vzero[1],vzero[2],
	      vzero[0],vzero[1],vzero[2],
	      izero); 
      else TheEVTtoWrite->AddMuonPositionInfo(MuonIds[ip],1,0.,0.,0.,0.,0.,0.,0.,0.);
    }
    //record stopping position//////////////
    if(!useANTARESformat)fprintf(outfile,"%.6e %.6e %.6e %.7e\n",
				 (stopPosition[ip])[0]/cm,(stopPosition[ip])[1]/cm,(stopPosition[ip])[2]/cm,(stopTime[ip])/ns); 
    else TheEVTtoWrite->AddMuonPositionInfo(MuonIds[ip],2,
					    (stopPosition[ip])[0]/m,(stopPosition[ip])[1]/m,(stopPosition[ip])[2]/m,
					    (stopTime[ip])/ns);
  }
  //write information of muon energies every 10 meters
#ifdef G4MYMUON_KEEPENERGY
  if(useANTARESformat){
    for(int ien=0 ; ien<EnergyAtPosition.size() ; ien++)EnergyAtPosition[ien]=EnergyAtPosition[ien]/GeV; //convert to GeV
    TheEVTtoWrite->AddMuonEnergyInfo(EnergyAtPosition);
  }
#endif
  //write to output file
  if(useANTARESformat)TheEVTtoWrite->WriteEvent();
#if defined(G4MYEM_PARAMETERIZATION) || defined(G4MYHA_PARAMETERIZATION) //newha
  long double cont;
  G4int TotalNumberOfCathods=MyStDetector->allCathods->GetNumberOfCathods();
  G4bool FineBin=false;
  G4int VertexSolidAngleBins=51;
  if(MyGenerator->ParamEnergy == 0.0){
    FineBin=true;
    VertexSolidAngleBins=71;
  }
  for(G4int icath=0; icath<TotalNumberOfCathods ; icath++){
    //    G4ThreeVector FromGeneToOM = MyStDetector->allCathods->GetPosition(icath) - MyGenerator->position;
    //G4double cosangle1=(MyGenerator->direction).dot(FromGeneToOM)/dist;

    G4double dist= MyStDetector->allCathods->GetCathodRadius(icath);
    G4int distbin=icath; //for the definition of distances look at the gdml geometry file
    G4int cang1bin;
    for(cang1bin=0 ; cang1bin <VertexSolidAngleBins ; cang1bin++){
      G4double weight=2*pi*dist*dist;
      if(!FineBin){
	if(cang1bin<15)weight*=0.1;
	else if(cang1bin<23)weight*=0.025;
	else if(cang1bin<43)weight*=0.005;
	else weight*=0.025;
      }
      else{
	if(cang1bin<15)weight*=0.1;
	else if(cang1bin<23)weight*=0.025;
	else if(cang1bin<28)weight*=0.005;
	else if(cang1bin<53)weight*=0.001;
	else if(cang1bin<63)weight*=0.005;
	else weight*=0.025;
      }
      // if(cosangle1<0.5){
      //   cang1bin=int((cosangle1+1.0)/0.1);} //bins 0-14
      // else if(cosangle1<0.7){
      //   cang1bin=int((cosangle1-0.5)/0.025)+15;} //bins 15-22
      // else if(cosangle1<0.8){
      //   cang1bin=int((cosangle1-0.7)/0.005)+23;} //bins 23-42
      // else if(cosangle1<1.0){
      //   cang1bin=int((cosangle1-0.8)/0.025)+43;} //bins 43-50
      // else cang1bin=50;                          //bin 50
      
      G4int ibin_d1=cang1bin+VertexSolidAngleBins*distbin;
      (*myCumNorma)[ibin_d1]+=weight;
      
      cont=(*myPhotonsNumber)[ibin_d1];//newbin
      (*myCumPhotons)[ibin_d1] += cont;//newbin
      (*myCumPhotonsRms)[ibin_d1] += cont*cont;//newbin
    }
  }
#endif
}

