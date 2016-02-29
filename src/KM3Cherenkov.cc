#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Poisson.hh"
#include "G4EmProcessSubType.hh"
#include "G4LossTableManager.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"
#include "KM3Cherenkov.hh"
#ifndef G4DISABLE_PARAMETRIZATION
#include "G4SDManager.hh"
#include "KM3SD.hh"
#ifdef G4TRACK_INFORMATION
#include "KM3TrackInformation.hh"
#endif
#endif

KM3Cherenkov::KM3Cherenkov(const G4String& processName, G4ProcessType type)
  : G4VProcess(processName, type)
{
  SetProcessSubType(fCerenkov);

  fTrackSecondariesFirst = false;
  fMaxBetaChange = 0.;
  fMaxPhotons = 0;

  thePhysicsTable = NULL;

  if (verboseLevel>0) {
    G4cout << GetProcessName() << " is created " << G4endl;
  }
  BuildThePhysicsTable();

  M_PI2=2*M_PI;
  MinMeanNumberOfPhotonsForParam=20.0;
#if !defined(G4MYEM_PARAMETERIZATION) && !defined(G4MYHA_PARAMETERIZATION)
#ifdef G4ENABLE_MIE
#ifndef G4DISABLE_PARAMETRIZATION
  poskeep = new std::vector<G4ThreeVector>;
  timekeep = new std::vector<G4double>;
  idprikeep = new std::vector<G4int>;
  depenekeep = new std::vector<G4double>;
  dirkeep = new std::vector<G4ThreeVector>;
#endif
#endif
#endif

#ifdef G4JUST_COUNT_PHOTONS
  Count_Photons=0.0;
  Posit_Photons_Mean=0.0;
  for(G4int i=0 ; i<3002 ; i++)Posit_Photons_Histo[i]=0.0;
#endif

}

KM3Cherenkov::~KM3Cherenkov() 
{
  if (thePhysicsTable != NULL) {
    thePhysicsTable->clearAndDestroy();
    delete thePhysicsTable;
  }
#if !defined(G4MYEM_PARAMETERIZATION) && !defined(G4MYHA_PARAMETERIZATION)
#ifdef G4ENABLE_MIE
#ifndef G4DISABLE_PARAMETRIZATION
  poskeep->clear();
  delete poskeep;
  timekeep->clear();
  delete timekeep;
  idprikeep->clear();
  delete idprikeep;
  depenekeep->clear();
  delete depenekeep;
  dirkeep->clear();
  delete dirkeep;
#endif
#endif
#endif

#ifndef G4DISABLE_PARAMETRIZATION
#ifdef G4ENABLE_MIE
  delete myFlux;
#endif
#endif

#ifdef G4JUST_COUNT_PHOTONS
  G4int ibin;
  long double cumul=0;
  for(ibin=0 ; ibin<3002 ; ibin++){
    cumul+=Posit_Photons_Histo[ibin];
    if(cumul>0.5*Count_Photons)break;
  }
  G4double Posit_Photons_Median=-10*m+0.5*cm+G4double(ibin-1)*cm;
  Posit_Photons_Mean/=Count_Photons;
  printf("Count_Photons %.20Le %.5Le %.5e\n",Count_Photons,Posit_Photons_Mean/m,Posit_Photons_Median/m);
#endif

}

//---------------------------------------------------------------------------------------------------------------------------

void KM3Cherenkov::SetDetector(KM3Detector* adet)
{
  MyStDetector = adet;
  MaxAbsDist=MyStDetector->MaxAbsDist;
#ifndef G4DISABLE_PARAMETRIZATION
#ifdef G4ENABLE_MIE
  myFlux =new KM3EMDirectFlux(MyStDetector->EMParametrization_FILE,MyStDetector->TotCathodArea);
#endif
#endif

}

#if !defined(G4ENABLE_MIE) || (defined(G4ENABLE_MIE) && !defined(G4DISABLE_PARAMETRIZATION)) //newmie
//the following function does exactly the oposite as the rotateUz. It rotates the vector x to a coordinate system that the p0 is pointing to the positive z-axis
//the rotation matrix is the transverse of the one at rotateUz. the vector p0 must be normalized !!!
void KM3Cherenkov::myrotate(G4ThreeVector &x, const G4ThreeVector& p0)
{
  G4double u1=p0.x();
  G4double u2=p0.y();
  G4double u3=p0.z();
  G4double up=u1*u1+u2*u2;
  G4double x0 = x.x();
  G4double x1 = x.y();
  G4double x2 = x.z();
  
  if (up>0) {
    up = sqrt(up);
    x.setX((u1*u3*x0 + u2*u3*x1)/up - up*x2);
    x.setY(  (-u2*x0   +  u1*x1)/up);
    x.setZ(    u1*x0   +  u2*x1   +   u3*x2);
  }
  else if (u3 < 0.) { x.setX(-x0); x.setZ(-x2); }      // phi=0  teta=pi
  else {};
}

G4int KM3Cherenkov::PhotonHitsaBenthos(G4double x1,G4double y1,G4double z1,G4double px,G4double py,G4double pz,G4double x0,G4double y0,G4double z0,G4double r, G4double dir1, G4double dir2, G4double dir3)
{
  G4double s2,s3,d,a1,xa1,ya1,za1;

  //Firstly lets check if the photon will cross the sphere!
  //Assuming (x-x0)^2+(y-y0)^2+(z-z0)^2=r^2 as the sphere equation and
  // (x-x1,y-y1,z-z1)=a(px,py,pz), where a is positive, the equation of the photon route.
  s2=(x1-x0)*px+(y1-y0)*py+(z1-z0)*pz;
  s3=(x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0)-r*r;
  d=s2*s2-s3;
  
  if(d>0)
  {		
    a1=-s2-sqrt(d); //these is the earliest solution of the equation a^2+2*s2*a+s3=0
    if(a1<0) {/*printf("The photon's is generated inside or it has already crossed the sphere!!\n");*/return(0);}
    else
    {
      return (1);  // we do not look if it hits up or down because if it falls on the back of the benthos it can be still detected
                   // may be we could try to lower the radius of the benthos at this point
                   // the code below is never executed
      za1=z1+a1*pz; //the z-coordinate of the earliest cross
      if((dir3==1) && (za1>z0)) return (1); //the photon firstly cross the active part of the photocathod (pmt looking up)
      else if((dir3==-1) && (za1<z0)) return (1); //the photon firstly cross the active part of the photocathod (pmt looking down)
      else if((dir3>-1) && (dir3<1))
	{
	  xa1=x1+a1*px;
	  ya1=y1+a1*py;
	  if( dir1*(xa1-x0)+dir2*(ya1-y0)+dir3*(za1-z0) >0 ) return (1);   //the photon firstly cross the active part of the photocathod (pmt looking sideways)
	  else return (0);
	}
      else return(0); //The photon firstly cross the dead part of the photocathod
    }
  }
  else return(0); //The photon does not cross the sphere neither at the future nor at the past. (d<=0).
}
//---------------------------------------------------------------------------------------------------------------------------
//This routine takes as input the vectors of the position and momentum of the photon and some strange number which might need
//to be wiped out, and ???? ###### for the hit benthos (from array HITBENTHOS) checks if the benthos are hit.
// ##################### UPDATE ########################################

G4int KM3Cherenkov::PhotonHitsAnyBenthos(G4ThreeVector r,G4ParticleMomentum p)
{
  //G4cout <<"x="  <<  r(0) << "y="  << r(1) << "z="  << r(2) << G4endl;
  //G4cout <<"px=" <<  p(0) << "py=" << p(1) << "pz=" << p(2) << G4endl;

  G4double xx,yy,zz,rr,dir1,dir2,dir3; 
  for (int i=0;i<(int)HITBENTHOS[0][0];i++)        //possibly hit benthos
    {
      if ((int)HITBENTHOS[i][9]==1)                       //check only the benthos for which phi (not visible here) is between the benthos minPhi and maxPhi
	{
	  xx=HITBENTHOS[i][3];           //get benthos x position
	  yy=HITBENTHOS[i][4];           //get benthos y position
	  zz=HITBENTHOS[i][5];           //get benthos z position
	  rr=HITBENTHOS[i][6];           //get benthos radius
	  if (PhotonHitsaBenthos(r(0),r(1),r(2),p(0),p(1),p(2),xx,yy,zz,rr,dir1,dir2,dir3)==1) return(1);
	}
    }
  return(0);
}

//---------------------------------------------------------------------------------------------------------------------------
//
//
G4int KM3Cherenkov::checkIfParticleCanEmitToShpere(G4ThreeVector center,G4double r, const G4double minCos, const G4double maxCos, G4double &minPhi, G4double &maxPhi, G4int icare)
{
  G4double maxCosTheta = -2;
  G4double minCosTheta =  2;
  G4double x,y,R,R2,RRpre,RRpost,r2,zcpre,zcpost;
  G4double cosphi,costheta,sinphi,sintheta;

  //in the following x,y are the coordinates of the center of the sphere so the sphere is at the xy-plane and the particle at the z-axis
  //zc is the z-axis position of the particle (on the z-axis)
  //the shift for this is (parent(0),parent(1),center(2)) for the prestep (parent) and poststep (parent1)

  // In the following we consider 2 cases. The one is when the trajectory of the particle penetrates the sphere and the other when it does not
  // In the second case the maximum angle is the angle it sees the front at the sphere from the poststep
  // the minimum angle is the angle it sees the back of the sphere from the prestep.
  // In the first case the minimum angle is 0 and the maximum is the same as in the second case if the particle is below the sphere.
  // If the particle is above the sphere the minimum angle is as in case 1, while the maximum is pi
  // keep in mind that max and min refer to the maximum and minimum angles (reverse for cosinus)
  
  zcpre = parent(2) - center(2);
  if(zcpre>r)return 0;         //the particle is in a position such that sees the sphere with angle greater then 90 degrees (not possible to emmit detectable photons)
  zcpost = parent1(2) - center(2);
  x = center(0) - parent(0);
  y = center(1) - parent(1);
  R2  = x*x + y*y;           //the distance^2 of the sphere's center from the center of the frame
  r2 = r*r;
  
  if(R2>r2)              //the second case
    {
      RRpre=R2+zcpre*zcpre;  //the distance^2 of the prestep from the center of the sphere
      R=sqrt(R2);
      //the following is cos(phi-theta) were phi=the angle the particle sees the center of the sphere at the prestep
      //                               and theta=the half-angle the particles sees the sphere (half opening)
      minCosTheta=(-zcpre*sqrt(RRpre-r2)+r*R)/RRpre;

      RRpost=R2+zcpost*zcpost;   //the distance^2 of the poststep from the center of the sphere
      //the following is cos(phi+theta) were phi=the angle the particle sees the center of the sphere at the poststep
      //                               and theta=the half-angle the particles sees the sphere (half opening)
      maxCosTheta=(-zcpost*sqrt(RRpost-r2)-r*R)/RRpost;

    }
  else                    //the first case
    {
      if(zcpre>0){        //the particle is on the upper half of the sphere
	RRpre=R2+zcpre*zcpre;
	if(RRpre<=r2){    //the particle transverse the sphere
	  if(icare==1){minPhi=0;maxPhi=M_PI2;}    //in this case if we speak of benthos the phi must contain all angles
	  return 1;
	}
	else{
	  R=sqrt(R2);
	  minCosTheta=(-zcpre*sqrt(RRpre-r2)+r*R)/RRpre;
	  maxCosTheta=-1;	  
	}
      }
      else if(zcpost<0){  //the particle is on the lower half of the sphere
	RRpost=R2+zcpost*zcpost;
	if(RRpost<=r2){    //the particle transverse the sphere
	  if(icare==1){minPhi=0;maxPhi=M_PI2;}    //in this case if we speak of benthos the phi must contain all angles
	  return 1;
	}
	else{
	  R=sqrt(R2);
	  maxCosTheta=(-zcpost*sqrt(RRpost-r2)-r*R)/RRpost;
	  minCosTheta=1;
	}
      }
      else{                                      //the particle is transversing the xy-plane either inside or tangentially to the sphere
        if(icare==1){minPhi=0;maxPhi=M_PI2;}    //in this case if we speak of benthos the phi must contain all angles
	return 1;
      }
    }
  
  if((minCos<maxCosTheta) || (maxCos>minCosTheta))    //
    {return (0);}
  else
    {
      if(icare==0) return(1);   //the icare is for the calculation of the phi values
      globalMaxCos = maxCosTheta;
      globalMinCos = minCosTheta;
      if(r>=R){
	minPhi = 0;
	maxPhi = M_PI2;}
      else{
	G4double phiGamma = atan2(y,x); 
	if (phiGamma<0) phiGamma = phiGamma + M_PI2;
	G4double dphiGamma=asin(r/R);
	minPhi = phiGamma - dphiGamma;
	maxPhi = phiGamma + dphiGamma;}
      return (1);
    }
}

//---------------------------------------------------------------------------------------------------------------------------

G4int KM3Cherenkov::mycheckParticleOneStar(const G4ThreeVector &p0, const G4ThreeVector &x0, const G4ThreeVector &xx0, const G4double &minCos, const G4double &maxCos)
{
  // first assume that the sphere contains the whole of the detector
  parent = x0 ;
  parent1 = xx0 ;

  myrotate(parent,p0);
  myrotate(parent1,p0);
  
  icountHitBenthos=0;
  Spheres* mySphere=MyStDetector->allSpheres;
  myIterativeCheck(mySphere,p0,minCos,maxCos);
  HITBENTHOS[0][0] = icountHitBenthos;
  return icountHitBenthos;
}

void KM3Cherenkov::myIterativeCheck(Spheres* mySphere,const G4ThreeVector &p0,const G4double &minCos,const G4double &maxCos)
{
  G4double minPhi,maxPhi;  
  G4ThreeVector centerkeep;
  G4int icare;
  G4double MaxDistFromCenter2;

  G4ThreeVector center=mySphere->center;
  G4double r=mySphere->radius;
  if(mySphere->allnext->size() == 0){centerkeep=center; icare =1; }  //end of the line
  else icare = 0;
  myrotate(center,p0);
  MaxDistFromCenter2=(r+MaxAbsDist)*(r+MaxAbsDist);
  if( MaxDistFromCenter2 < (center-parent).mag2() ) return;   
  if(checkIfParticleCanEmitToShpere(center,r,minCos,maxCos,minPhi,maxPhi,icare)==0) return;  //here we check only the angles of Cherenkov emission
  if(icare==1){
    HITBENTHOS[icountHitBenthos][1] = minPhi;
    HITBENTHOS[icountHitBenthos][2] = maxPhi;
    HITBENTHOS[icountHitBenthos][3] = centerkeep(0);
    HITBENTHOS[icountHitBenthos][4] = centerkeep(1);
    HITBENTHOS[icountHitBenthos][5] = centerkeep(2);
    HITBENTHOS[icountHitBenthos][6] = r;
    HITBENTHOS[icountHitBenthos][7] = globalMaxCos;
    HITBENTHOS[icountHitBenthos][8] = globalMinCos;
    HITBENTHOS[icountHitBenthos][9] = 0.0;
    icountHitBenthos++;
    return;
  }
  for(size_t isp=0 ; isp<mySphere->allnext->size() ; isp++){
    myIterativeCheck( (*(mySphere->allnext))[isp],p0,minCos,maxCos);
  }
}

//---------------------------------------------------------------------------------------------------------------------------
//Given a phi angle this routine checks (for the possibly hit benthos as given by HITBENTHOS[0][0]) whether this angle falls
//within the range defined by minPhi (HITBENTHOS[i][1]) and maxPhi (HITBENTHOS[i][2]). In such a case the routine returns 1,
//otherwise a zero.


G4int KM3Cherenkov::checkPhi(const G4double& aphi)
{
  G4int flag=0;
  for (int i=0;i<(int)HITBENTHOS[0][0];i++)
    {
      HITBENTHOS[i][9]=0.0;
      if ((HITBENTHOS[i][1]<=aphi) && (aphi<=HITBENTHOS[i][2]))
	{
	  HITBENTHOS[i][9]=1;
	  flag=1;
	}
      else if (HITBENTHOS[i][1]<0)
	{
	  if ((HITBENTHOS[i][1]+M_PI2<=aphi) || (aphi<=HITBENTHOS[i][2]))
	    {
	      HITBENTHOS[i][9]=1;
	      flag=1;
	    }
	}
      else if  (HITBENTHOS[i][2]>M_PI2)
	{
	  if ((HITBENTHOS[i][1]<=aphi) || (aphi<=HITBENTHOS[i][2]-M_PI2))
	    {
	      HITBENTHOS[i][9]=1;
	      flag=1;
	    }
	}
    }
  return(flag);
}

//---------------------------------------------------------------------------------------------------------------------------

#endif

//---------------------------------------------------------------------------------------------------------------------------

G4VParticleChange* KM3Cherenkov::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)

// This routine is called for each tracking Step of a charged particle
// in a radiator. A Poisson-distributed number of photons is generated
// according to the Cerenkov formula, distributed evenly along the track
// segment and uniformly azimuth w.r.t. the particle direction. The 
// parameters are then transformed into the Master Reference System, and 
// they are added to the particle change. 

{
  //G4cout << " IN CHERENKOV";
  //G4cout << " Particle id " << aTrack.GetTrackID() << G4endl;
  //G4cout << " Parent id " << aTrack.GetParentID();
  //G4cout << " Particle name " << aTrack.GetDefinition()->GetParticleName();
  //G4cout << " Kinetic Energy " << aTrack.GetKineticEnergy() << " " << aTrack.GetDynamicParticle()->GetKineticEnergy();
  //G4cout << " Total Energy " << aTrack.GetTotalEnergy();
  //G4cout << " Creator Process " << aTrack.GetCreatorProcess()->GetProcessName();
  //G4cout << " Track Length " << aTrack.GetTrackLength();
  //G4cout << " STEP Length from step " << aStep.GetStepLength();
  //G4cout << " delta position phi " << aStep.GetDeltaPosition().phi()*deg;
  //G4cout << " delta position theta " << aStep.GetDeltaPosition().theta()*deg;
  //G4cout << " delta position magnitude " << aStep.GetDeltaPosition().mag()<<G4endl;

  ///at first initialize the pointers to Q_E, glass and gell transparencies////
  static G4MaterialPropertyVector* QECathod = NULL;
#ifdef G4MY_TRANSPARENCIES
  static G4MaterialPropertyVector* AbsBenth = NULL;
  static G4MaterialPropertyVector* AbsGell = NULL;
#endif
  if(QECathod == NULL){
    const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
    for (size_t J=0 ; J<theMaterialTable->size() ; J++) {
      if ((*theMaterialTable)[J]->GetName() == G4String("Cathod")){
	G4MaterialPropertiesTable* aMaterialPropertiesTable = (*theMaterialTable)[J]->GetMaterialPropertiesTable();
	QECathod=aMaterialPropertiesTable->GetProperty("Q_EFF");
      }
#ifdef G4MY_TRANSPARENCIES
      else if ((*theMaterialTable)[J]->GetName() == G4String("Glass")){
	G4MaterialPropertiesTable* aMaterialPropertiesTable = (*theMaterialTable)[J]->GetMaterialPropertiesTable();
	AbsBenth=aMaterialPropertiesTable->GetProperty("ABSLENGTH");
      }
      else if((*theMaterialTable)[J]->GetName() == G4String("Gell")){
	G4MaterialPropertiesTable* aMaterialPropertiesTable = (*theMaterialTable)[J]->GetMaterialPropertiesTable();
	AbsGell=aMaterialPropertiesTable->GetProperty("ABSLENGTH");
      }
#endif
    }
  }
  ///end of initialization//////////////////////////////////////////////////////

  //////////////////////////////////////////////////////
  // Should we ensure that the material is dispersive?
  //////////////////////////////////////////////////////
  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4Material* aMaterial = aTrack.GetMaterial();

  G4MaterialPropertiesTable* aMaterialPropertiesTable =
    aMaterial->GetMaterialPropertiesTable();

  if (!aMaterialPropertiesTable)
    return pParticleChange;

  G4MaterialPropertyVector* Rindex = 
    aMaterialPropertiesTable->GetProperty("RINDEX"); 

  if (!Rindex) 
    return pParticleChange;

#ifdef G4MYHAMUONS_PARAMETERIZATION
  aParticleChange.SetNumberOfSecondaries(0);
  return pParticleChange;
#endif

  //check that the particle is inside the active volume of the detector
  static G4double detectorMaxRho2=MyStDetector->detectorMaxRho * MyStDetector->detectorMaxRho;
  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4ThreeVector x0 = pPreStepPoint->GetPosition();
  //  G4cout <<"prepoint "<< x0[0] <<" "<< x0[1] <<" "<< x0[2] <<G4endl;
  G4ThreeVector distanceV=x0 - MyStDetector->detectorCenter;
  G4double distanceRho2=distanceV[0]*distanceV[0] + distanceV[1]*distanceV[1];
  if( (distanceRho2>detectorMaxRho2) || (x0[2]<MyStDetector->bottomPosition) || (x0[2]>MyStDetector->detectorMaxz)) {
    // return unchanged particle and no secondaries
    aParticleChange.SetNumberOfSecondaries(0);
    return pParticleChange;
  }

  //step length
  G4double step_length = aStep.GetStepLength();
  //  G4cout << "step_length "<<aParticle->GetDefinition()->GetParticleName() <<" "<< step_length<<G4endl;
  // particle charge
  const G4double charge = aParticle->GetDefinition()->GetPDGCharge();
  // particle beta
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
  const G4double beta = (pPreStepPoint ->GetBeta() +
			 pPostStepPoint->GetBeta())/2.;
  G4double MeanNumPhotons = GetAverageNumberOfPhotons(charge,beta,aMaterial,Rindex);
  MeanNumPhotons *= step_length*MyStDetector->Quantum_Efficiency;
  G4double BetaInverse = 1.0/beta;
  G4ThreeVector p0 = aStep.GetDeltaPosition().unit();

  G4bool EmittedAsScattered=false; //newmie

#ifndef G4DISABLE_PARAMETRIZATION
#ifdef G4ENABLE_MIE
  //direct photons global parametrization section only for primary muons
  if( ((aTrack.GetDefinition()->GetParticleName() == "mu-") || (aTrack.GetDefinition()->GetParticleName() == "mu+")) && aTrack.GetParentID()==0 ){
    if(BetaInverse < 1.01){ //this accounts for v>0.99c
      //      G4cout <<"StepLength "<<step_length/cm<<G4endl;
      poskeep->push_back(x0);
      timekeep->push_back(pPreStepPoint->GetGlobalTime());
      idprikeep->push_back(aTrack.GetTrackID());
      depenekeep->push_back(MeanNumPhotons);
      //      G4cout<<"PostStepDoIt "<<x0[0]/m<<" "<<x0[1]/m<<" "<<x0[2]/m<<" "<<pPreStepPoint->GetGlobalTime()<<" "<<MeanNumPhotons<<G4endl;
      dirkeep->push_back(p0);
      // return unchanged particle and no secondaries
      //newmie aParticleChange.SetNumberOfSecondaries(0);
      //newmie return pParticleChange;
      EmittedAsScattered=true; //newmie
    }
  }
  ///////////////////////////////////////////////
#endif
#endif

  if (MeanNumPhotons <= 0.0) {
    // return unchanged particle and no secondaries
    aParticleChange.SetNumberOfSecondaries(0);
    return pParticleChange;
  }

  G4int NumPhotons = (G4int) CLHEP::RandPoisson::shoot(MeanNumPhotons);
  if (NumPhotons <= 0) {
    // return unchanged particle and no secondaries  
    aParticleChange.SetNumberOfSecondaries(0);
    return pParticleChange;
  }

#ifdef G4JUST_COUNT_PHOTONS
  Count_Photons+=(long double)NumPhotons;
  G4double pos_z_projection=0.5*( (aStep.GetPreStepPoint()->GetPosition())[2]+
				  (aStep.GetPostStepPoint()->GetPosition())[2]);
  Posit_Photons_Mean+=((long double)NumPhotons)*((long double)pos_z_projection);
  G4int ibin=1001+int(floor(pos_z_projection/cm));
  if(ibin<0)ibin=0;
  if(ibin>3001)ibin=3001;
  Posit_Photons_Histo[ibin]+=(long double)NumPhotons;
  aParticleChange.SetNumberOfSecondaries(0);
  return pParticleChange;
#endif

#ifndef G4DISABLE_PARAMETRIZATION
#ifdef G4ENABLE_MIE
  //here is the parametrization step oriented only for a minimum number of secondaries
  //initialization for the group speed at max qe
  static G4double thespeedmaxQE=-1.0;
  static KM3SD * aMySD = NULL;
  if (thespeedmaxQE<0){
    G4Material* aMaterial = G4Material::GetMaterial("Cathod");
    G4double MaxQE=-1;
    G4double PhEneAtMaxQE;
    G4MaterialPropertyVector* aPropertyVector = aMaterial->GetMaterialPropertiesTable()->GetProperty("Q_EFF");
    for (size_t i = 0; i < aPropertyVector->GetVectorLength() ; i++){
      G4double ThisQE=(*aPropertyVector)[i];
      G4double ThisPhEne=aPropertyVector->Energy(i);
      if(ThisQE > MaxQE){MaxQE=ThisQE;PhEneAtMaxQE=ThisPhEne;}
    }
    aMaterial = G4Material::GetMaterial("Water");
    G4MaterialPropertyVector* GroupVel = aMaterial->GetMaterialPropertiesTable()->GetProperty("GROUPVEL");
    thespeedmaxQE=GroupVel->Value(PhEneAtMaxQE); //coresponds to the maximum qe each time. This is the right one
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    aMySD = (KM3SD*)SDman->FindSensitiveDetector(G4String("mydetector1/MySD"),true);
  }
  //end of initialization
  if( (MeanNumPhotons > MinMeanNumberOfPhotonsForParam) && BetaInverse<1.01 && !EmittedAsScattered){ //this accounts for v>0.99c //newmie
    const G4VProcess * theProcess = aTrack.GetCreatorProcess();
    G4int originalTrackCreatorProcess;
    G4int originalParentID;
#ifdef G4TRACK_INFORMATION
    if(theProcess != NULL){
      KM3TrackInformation* info= (KM3TrackInformation*)(aTrack.GetUserInformation());
      originalParentID=info->GetOriginalParentID();
      G4String creator=info->GetOriginalTrackCreatorProcess();
      if(creator=="KM3Cherenkov")originalTrackCreatorProcess=0;
      else if(creator=="muPairProd")originalTrackCreatorProcess= 1;
      else if(creator=="muIoni")originalTrackCreatorProcess= 2;
      else if(creator=="muBrems")originalTrackCreatorProcess= 3;
      else if(creator=="muonNuclear")originalTrackCreatorProcess= 4;
      else if(creator=="Decay")originalTrackCreatorProcess= 8;
      else if(creator=="muMinusCaptureAtRest")originalTrackCreatorProcess= 9;
      else originalTrackCreatorProcess= 5;
    }
    else{ //it is an initial particle
      originalParentID=aTrack.GetTrackID();
      originalTrackCreatorProcess=0;  //it is cherenkov
    }
#else
    originalParentID=1;
    originalTrackCreatorProcess=0;
#endif
    G4int originalInfo = (originalParentID-1)*10 + originalTrackCreatorProcess ;
	  
    G4ThreeVector MiddleOfStep=0.5*(pPreStepPoint->GetPosition()+pPostStepPoint->GetPosition());
    G4double MiddleOfTime=0.5*(pPreStepPoint->GetGlobalTime()+pPostStepPoint->GetGlobalTime());
    //from here is the actual parametrization////////////////////////////////////////////////////////
    static G4double MaxAbsDist2=MyStDetector->MaxAbsDist * MyStDetector->MaxAbsDist;
    static G4int TotalNumberOfTowers=MyStDetector->allTowers->size();
    for(int it=0;it<TotalNumberOfTowers;it++){
      G4double dx=(*(MyStDetector->allTowers))[it]->position[0]-MiddleOfStep[0];
      G4double dy=(*(MyStDetector->allTowers))[it]->position[1]-MiddleOfStep[1];
      G4double distancetower2=dx*dx+dy*dy;
      if(distancetower2<MaxAbsDist2){
	G4int TotalNumberOfOMs=(*(MyStDetector->allTowers))[it]->BenthosIDs->size();
	for(int iot=0;iot<TotalNumberOfOMs;iot++){
	  G4int io=(*(*(MyStDetector->allTowers))[it]->BenthosIDs)[iot];
	  G4ThreeVector FromGeneToOM=(*MyStDetector->allOMs)[io]->position - MiddleOfStep;
	  G4double distancein=FromGeneToOM.mag2();
	  if(distancein<MaxAbsDist2){
	    distancein=sqrt(distancein);
	    FromGeneToOM /= distancein;
	    G4double anglein=p0.dot(FromGeneToOM);
	    myFlux->FindBins(MeanNumPhotons,distancein,anglein);  //here change
	    G4int NumberOfSamples=myFlux->GetNumberOfSamples();//here change
	    G4int icstart,icstop;
	    G4double theFastTime;
	    G4ThreeVector x,y,z;
	    if(NumberOfSamples>0){
	      icstart=(*(*MyStDetector->allOMs)[io]->CathodsIDs)[0];
	      icstop=1+(*(*MyStDetector->allOMs)[io]->CathodsIDs)[(*MyStDetector->allOMs)[io]->CathodsIDs->size() - 1];
	      theFastTime = distancein/thespeedmaxQE +MiddleOfTime;
	      z=FromGeneToOM;
	      y=p0.cross(z)/sqrt(1.0-anglein*anglein);
	      x=y.cross(z);
	    }
	    for(G4int isa=0 ; isa<NumberOfSamples ; isa++){
	      onePE aPE=myFlux->GetSamplePoint();//here change
	      //	G4cout << "OutFromParam "<<distancein<<" "<<anglein<<" "<<aPE.costh<<" "<<aPE.phi<<" "<<aPE.time<<G4endl;  //tempo
	      G4double costh=aPE.costh;//here change
	      G4double sinth=sqrt(1.0-costh*costh);
	      G4double cosphi=cos(aPE.phi);//here change
	      G4double sinphi=sin(aPE.phi);//here change
	      //short	      G4ThreeVector photonDirection=-(sinth*(cosphi*x+sinphi*y)+costh*z);
	      G4ThreeVector photonDirection=(sinth*(cosphi*x+sinphi*y)+costh*z);
	      //short	      G4double angleThetaDirection=photonDirection.theta();
	      //short	      G4double anglePhiDirection=photonDirection.phi();
	      //short	      angleThetaDirection *= 180./M_PI;
	      //short	      anglePhiDirection *= 180./M_PI;
	      //short	      if(anglePhiDirection < 0.0)anglePhiDirection += 360.0;
	      //short	      G4int angleDirection=(G4int)(nearbyint(angleThetaDirection)*1000.0 + nearbyint(anglePhiDirection));
	      G4int ic=G4int(icstart+(icstop-icstart)*G4UniformRand());
	      //short	      aMySD->InsertExternalHit(ic,theFastTime+aPE.time,originalInfo,angleDirection,-998);//here change
	      aMySD->InsertExternalHit(ic,(*MyStDetector->allOMs)[io]->position,
				       theFastTime+aPE.time,originalInfo,photonDirection);//here change
	    }//for(G4int isa=0 ; isa<NumberOfSamples ; isa++){
	  }//if(distancein<MyStDetector->MaxAbsDist){
	}//for(int io=0;io<TotalNumberOfOMs;io++){
      }//if(distancetower2<MaxAbsDist2)
    }//for(int it=0;it<TotalNumberOfTowers;it++)
    ////////////////////////////////////////////////////////////////////////////////////////////////
    //at the end produce nothing for tracking
    //newmie aParticleChange.SetNumberOfSecondaries(0);
    //newmie return pParticleChange;
    EmittedAsScattered=true; //newmie
  } //if( (MeanNumPhotons > MinMeanNumberOfPhotonsForParam) && BetaInverse<1.01){
  //////////////////////////////////////////////////////////////////////////////
#endif
#endif
  G4double nMax = Rindex->GetMaxValue();
  G4double maxCos = BetaInverse / nMax; 

#if !defined(G4ENABLE_MIE) || (defined(G4ENABLE_MIE) && !defined(G4DISABLE_PARAMETRIZATION)) //newmie
  //here is the classic iterative check for no scattering///////////////////////////////////////////
  ////////Check the if the charged particle can emit to any benthos
  G4ThreeVector xx0=pPostStepPoint->GetPosition();
  G4double nMin = Rindex->GetMinValue();
  G4double minCos = BetaInverse / nMin; 
#if !defined(G4ENABLE_MIE)  //newmie
  if(MeanNumPhotons > MinMeanNumberOfPhotonsForParam){
#else
  if(MeanNumPhotons > MinMeanNumberOfPhotonsForParam && EmittedAsScattered){ //newmie
#endif
    G4int numben=0;
    numben=mycheckParticleOneStar(p0,x0,xx0,minCos,maxCos);
    if (numben==0) {
      // return unchanged particle and no secondaries  
      aParticleChange.SetNumberOfSecondaries(0);
      return pParticleChange;
    }
  }
  ////end of check
  ///////////////////////////////////////////////////////////////////////////////////////////////////
#endif
	
  aParticleChange.SetNumberOfSecondaries(NumPhotons); 

  if (fTrackSecondariesFirst)
    {
      if (aTrack.GetTrackStatus() == fAlive )
	aParticleChange.ProposeTrackStatus(fSuspend);
    }

  G4double Pmin = Rindex->GetMinLowEdgeEnergy();
  G4double Pmax = Rindex->GetMaxLowEdgeEnergy();
  G4double dp = Pmax - Pmin;
  G4double maxSin2 = (1.0 - maxCos) * (1.0 + maxCos);
  const G4double beta1 = pPreStepPoint ->GetBeta();
  const G4double beta2 = pPostStepPoint->GetBeta();
  G4double MeanNumberOfPhotons1 =
    GetAverageNumberOfPhotons(charge,beta1,aMaterial,Rindex);
  G4double MeanNumberOfPhotons2 =
    GetAverageNumberOfPhotons(charge,beta2,aMaterial,Rindex);
  G4double t0 = pPreStepPoint->GetGlobalTime();
  //  NumPhotons=0; //lookout
  for (G4int i = 0; i < NumPhotons; i++){
    
    G4double rand;
    G4double sampledEnergy, sampledRI; 
    G4double cosTheta, sin2Theta;
    
    //sample a phi
    rand=G4UniformRand();      
    G4double phi=M_PI2*rand;
#if !defined(G4ENABLE_MIE) || (defined(G4ENABLE_MIE) && !defined(G4DISABLE_PARAMETRIZATION)) //newmie
    //here check if this phi can hit any benthos
    G4int chch;
#if !defined(G4ENABLE_MIE) //newmie
    if(MeanNumPhotons > MinMeanNumberOfPhotonsForParam)chch=checkPhi(phi);
#else
    if(MeanNumPhotons > MinMeanNumberOfPhotonsForParam && EmittedAsScattered)chch=checkPhi(phi); //newmie
#endif
    else chch=1;
    if (chch==1){             
#endif
      G4double sinPhi = sin(phi);
      G4double cosPhi = cos(phi);
      
      // Determine photon energy
      // sample an energy
      
      do {
	rand = G4UniformRand();     
	sampledEnergy = Pmin + rand * dp; 
	sampledRI = Rindex->Value(sampledEnergy);
	cosTheta = BetaInverse / sampledRI;  
	
	sin2Theta = (1.0 - cosTheta)*(1.0 + cosTheta);
	rand = G4UniformRand();	
	
      } while (rand*maxSin2 > sin2Theta);
      
      G4double qeProb=QECathod->Value(sampledEnergy);
#ifdef G4MY_TRANSPARENCIES
      qeProb *= exp(-15.0/AbsBenth->Value(sampledEnergy) - 10.0/AbsGell->Value(sampledEnergy));
#endif
      if (G4UniformRand()<qeProb){
	// calculate x,y, and z components of photon momentum 
	// (in coord system with primary particle direction 
	//  aligned with the z axis)
	
	G4double sinTheta = sqrt(sin2Theta); 
	G4double px = sinTheta*cosPhi;
	G4double py = sinTheta*sinPhi;
	G4double pz = cosTheta;
	
	// Create photon momentum direction vector 
	// The momentum direction is still with respect
	// to the coordinate system where the primary
	// particle direction is aligned with the z axis  
	
	G4ParticleMomentum photonMomentum(px, py, pz);
	
	// Rotate momentum direction back to global reference
	// system 
	
	photonMomentum.rotateUz(p0);
	
	// Determine polarization of new photon 
	
	G4double sx = cosTheta*cosPhi;
	G4double sy = cosTheta*sinPhi; 
	G4double sz = -sinTheta;
	
	G4ThreeVector photonPolarization(sx, sy, sz);
	
	// Rotate back to original coord system 
	
	photonPolarization.rotateUz(p0);
	
	// Generate new G4Track object and a new photon
	//first find generation position and time
	
	G4double delta, NumberOfPhotons, N;
	do {
	  rand = G4UniformRand();
	  delta = rand * step_length;
	  NumberOfPhotons = MeanNumberOfPhotons1 - delta *
	    (MeanNumberOfPhotons1-MeanNumberOfPhotons2)/step_length;
	  N = G4UniformRand() *
	    std::max(MeanNumberOfPhotons1,MeanNumberOfPhotons2);
	} while (N > NumberOfPhotons);

	G4double deltaTime = delta /
	  ((pPreStepPoint->GetVelocity()+
	    pPostStepPoint->GetVelocity())/2.);
	
	G4ThreeVector aSecondaryPosition =
	  x0 + rand * aStep.GetDeltaPosition();
	
	G4double aSecondaryTime = t0 + deltaTime;
	
#if !defined(G4ENABLE_MIE) || (defined(G4ENABLE_MIE) && !defined(G4DISABLE_PARAMETRIZATION)) //newmie
	//here see if this photon can hit any benthos
	G4int  Bhit1;
#if !defined(G4ENABLE_MIE) //newmie
	if(MeanNumPhotons > MinMeanNumberOfPhotonsForParam)Bhit1=PhotonHitsAnyBenthos(aSecondaryPosition,photonMomentum);
#else
	if(MeanNumPhotons > MinMeanNumberOfPhotonsForParam && EmittedAsScattered)Bhit1=PhotonHitsAnyBenthos(aSecondaryPosition,photonMomentum); //newmie
#endif
	else Bhit1=1;
	if ((Bhit1>=1)){
#endif
	  
	  G4DynamicParticle* aCerenkovPhoton =
	    new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), 
				  photonMomentum);
	  aCerenkovPhoton->SetPolarization
	    (photonPolarization.x(),
	     photonPolarization.y(),
	     photonPolarization.z());
	  
	  aCerenkovPhoton->SetKineticEnergy(sampledEnergy);
	  
	  
	  // Generate the track
	  G4Track* aSecondaryTrack = 
	    new G4Track(aCerenkovPhoton,aSecondaryTime,aSecondaryPosition);
	  
	  aSecondaryTrack->SetTouchableHandle(aStep.GetPreStepPoint()->GetTouchableHandle());
	  aSecondaryTrack->SetParentID(aTrack.GetTrackID());
	  //newmie
#if defined(G4TRACK_INFORMATION) && !defined(G4DISABLE_PARAMETRIZATION)
	  if(EmittedAsScattered){
	    KM3TrackInformation * anInfo = new KM3TrackInformation();
	    aSecondaryTrack->SetUserInformation(anInfo);
	  }
#endif
	  //newmie
	  aParticleChange.AddSecondary(aSecondaryTrack);
#if !defined(G4ENABLE_MIE) || (defined(G4ENABLE_MIE) && !defined(G4DISABLE_PARAMETRIZATION)) //newmie
	}
#endif
      }// if (G4UniformRand()<qeProb){
#if !defined(G4ENABLE_MIE) || (defined(G4ENABLE_MIE) && !defined(G4DISABLE_PARAMETRIZATION)) //newmie
    }//if check phi
#endif
  }//for each photon
	

  if (verboseLevel>0) {
    G4cout << "\n Exiting from KM3Cherenkov::DoIt -- NumberOfSecondaries = " 
	   << aParticleChange.GetNumberOfSecondaries() << G4endl;
  }
	
  return pParticleChange;
}

#if !defined(G4MYEM_PARAMETERIZATION) && !defined(G4MYHA_PARAMETERIZATION)
#ifdef G4ENABLE_MIE
#ifndef G4DISABLE_PARAMETRIZATION
void KM3Cherenkov::CreateDirectPhotons()
{
  //here is the parametrization of direct photons from muons
  //initialization for the group speed at max qe
  static G4double thespeedmaxQE=-1.0;
  static KM3SD * aMySD = NULL;
  if (thespeedmaxQE<0){
    G4Material* aMaterial = G4Material::GetMaterial("Cathod");
    G4double MaxQE=-1;
    G4double PhEneAtMaxQE;
    G4MaterialPropertyVector* aPropertyVector = aMaterial->GetMaterialPropertiesTable()->GetProperty("Q_EFF");
    for (size_t i = 0; i < aPropertyVector->GetVectorLength() ; i++){
      G4double ThisQE=(*aPropertyVector)[i];
      G4double ThisPhEne=aPropertyVector->Energy(i);
      if(ThisQE > MaxQE){MaxQE=ThisQE;PhEneAtMaxQE=ThisPhEne;}
    }
    aMaterial = G4Material::GetMaterial("Water");
    G4MaterialPropertyVector* GroupVel = aMaterial->GetMaterialPropertiesTable()->GetProperty("GROUPVEL");
    thespeedmaxQE=GroupVel->Value(PhEneAtMaxQE); //coresponds to the maximum qe each time. This is the right one
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    aMySD = (KM3SD*)SDman->FindSensitiveDetector(G4String("mydetector1/MySD"),true);
  }
  //end of initialization
  static G4int originalTrackCreatorProcess=0; //it is always Direct Cherenkov
  static G4int TotalNumberOfTowers=MyStDetector->allTowers->size();
  static G4double MaxAbsDist2=MyStDetector->MaxAbsDist * MyStDetector->MaxAbsDist;
  static G4double distbin2=(50.0*cm)*(50.0*cm); //is the distance binning in the direct photons for gathering
  size_t arraysize=idprikeep->size();
  //  G4cout << "----------- size from direct ----------- "<<arraysize<<G4endl;
  size_t counter=0;
  while (counter<arraysize){
    G4int idpri=(*idprikeep)[counter];
    while( (counter<arraysize) && (idpri == (*idprikeep)[counter]) ){
      G4int originalInfo = ((*idprikeep)[counter] - 1)*10 + originalTrackCreatorProcess ; //to write in MySD
      G4ThreeVector pospri=(*poskeep)[counter];
      G4double timepri=(*timekeep)[counter];
      G4double depene=0.0;
      while( (counter<arraysize) && (idpri == (*idprikeep)[counter]) && ((pospri - (*poskeep)[counter]).mag2() < distbin2) ){
	depene += (*depenekeep)[counter];
	//	G4cout<<"CreateDirectPhotons1 "<<counter<<" "<<(*depenekeep)[counter]<<G4endl;
	counter++;
      }
      G4ThreeVector p0;
      G4ThreeVector thispos;
      G4double thistime;
      if( (counter == arraysize) || ((*idprikeep)[counter] != idpri) ){
	p0 = (*poskeep)[counter-1] - pospri;
	thispos=0.5*( (*poskeep)[counter-1] + pospri );
	//	G4cout<<"CreateDirectPhotons2 "<<counter-1<<" "<<(*poskeep)[counter-1]<<" "<<pospri<<G4endl;
	thistime=0.5*( (*timekeep)[counter-1] + timepri);
	//	G4cout<<"CreateDirectPhotons3 "<<counter-1<<" "<<(*timekeep)[counter-1]<<" "<<timepri<<G4endl; 
	G4double step=p0.mag();
	if(step == 0.0){
	  //before	  p0=pospri - (*poskeep)[counter-2];
	  p0 = (*dirkeep)[counter-1];  //this the case when one distance interval is composed by only one step at the end of the track
	}
      }
      else{
	p0 = (*poskeep)[counter] - pospri;
	thispos=0.5*( (*poskeep)[counter] + pospri );
	//	G4cout<<"CreateDirectPhotons4 "<<counter<<" "<<(*poskeep)[counter]<<" "<<pospri<<G4endl;
	thistime=0.5*( (*timekeep)[counter] + timepri);
	//	G4cout<<"CreateDirectPhotons5 "<<counter<<" "<<(*timekeep)[counter]<<" "<<timepri<<G4endl; 
      }
      //      G4cout<<"CreateDirectPhotons6 "<<thispos[0]/m<<" "<<thispos[1]/m<<" "<<thispos[2]/m<<" "<<thistime<<" "<<depene<<G4endl;
      p0 = p0.unit();
      for(int it=0;it<TotalNumberOfTowers;it++){
	G4double dx=(*(MyStDetector->allTowers))[it]->position[0]-thispos[0];
	G4double dy=(*(MyStDetector->allTowers))[it]->position[1]-thispos[1];
	G4double distancetower2=dx*dx+dy*dy;
	if(distancetower2<MaxAbsDist2){
	  G4int TotalNumberOfOMs=(*(MyStDetector->allTowers))[it]->BenthosIDs->size();
	  for(int iot=0;iot<TotalNumberOfOMs;iot++){
	    G4int io=(*(*(MyStDetector->allTowers))[it]->BenthosIDs)[iot];
	    G4ThreeVector FromGeneToOM=(*MyStDetector->allOMs)[io]->position - thispos;
	    G4double distancein=FromGeneToOM.mag2();
	    if(distancein<MaxAbsDist2){
	      distancein=sqrt(distancein);
	      FromGeneToOM /= distancein;
	      G4double anglein=p0.dot(FromGeneToOM);
	      myFlux->FindBins(depene,distancein,anglein);  //here change
	      G4int NumberOfSamples=myFlux->GetNumberOfSamples();//here change
	      G4int icstart,icstop;
	      G4double theFastTime;
	      G4ThreeVector x,y,z;
	      if(NumberOfSamples>0){
		icstart=(*(*MyStDetector->allOMs)[io]->CathodsIDs)[0];
		icstop=1+(*(*MyStDetector->allOMs)[io]->CathodsIDs)[(*MyStDetector->allOMs)[io]->CathodsIDs->size() - 1];
		theFastTime = distancein/thespeedmaxQE + thistime;
		z=FromGeneToOM;
		y=p0.cross(z)/sqrt(1.0-anglein*anglein);
		x=y.cross(z);
	      }
	      for(G4int isa=0 ; isa<NumberOfSamples ; isa++){
		onePE aPE=myFlux->GetSamplePoint();
		G4double costh=aPE.costh;//here change
		G4double sinth=sqrt(1.0-costh*costh);
		G4double cosphi=cos(aPE.phi);//here change
		G4double sinphi=sin(aPE.phi);//here change
		//		G4cout <<"InC "<<x.mag2()<<" "<<y.mag2()<<" "<<z.mag2()<<G4endl;
		//short		G4ThreeVector photonDirection=-(sinth*(cosphi*x+sinphi*y)+costh*z);
		G4ThreeVector photonDirection=(sinth*(cosphi*x+sinphi*y)+costh*z);
		//short		G4double angleThetaDirection=photonDirection.theta();
		//short		G4double anglePhiDirection=photonDirection.phi();
		//short		angleThetaDirection *= 180./M_PI;
		//short		anglePhiDirection *= 180./M_PI;
		//short		if(anglePhiDirection < 0.0)anglePhiDirection += 360.0;
		//short		G4int angleDirection=(G4int)(nearbyint(angleThetaDirection)*1000.0 + nearbyint(anglePhiDirection));
		G4int ic=G4int(icstart+(icstop-icstart)*G4UniformRand());
		//short		aMySD->InsertExternalHit(ic,theFastTime+aPE.time,originalInfo,angleDirection,-996);
		aMySD->InsertExternalHit(ic,(*MyStDetector->allOMs)[io]->position,
					 theFastTime+aPE.time,originalInfo,photonDirection);
	      }//for(G4int isa=0 ; isa<NumberOfSamples ; isa++){
	    }//if(distancein<MyStDetector->MaxAbsDist){
	  }//for(int io=0;io<TotalNumberOfOMs;io++){
	}//if(distancetower2<MaxAbsDist2)
      }//for(int it=0;it<TotalNumberOfTowers;it++)
    }//while( idpri == (*idprikeep)[counter] ){
  }//while (counter<arraysize){
  poskeep->clear();
  timekeep->clear();
  idprikeep->clear();
  depenekeep->clear();
  dirkeep->clear();
}
#endif
#endif
#endif


//---------------------------------------------------------------------------------------------------------------------------

void KM3Cherenkov::BuildThePhysicsTable()
{
  if (thePhysicsTable) return;

  const G4MaterialTable* theMaterialTable=
    G4Material::GetMaterialTable();
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();

  // create new physics table
	
  thePhysicsTable = new G4PhysicsTable(numOfMaterials);

  // loop for materials

  for (G4int i=0 ; i < numOfMaterials; i++)
    {
      G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector =
	new G4PhysicsOrderedFreeVector();

      // Retrieve vector of refraction indices for the material
      // from the material's optical properties table 

      G4Material* aMaterial = (*theMaterialTable)[i];

      G4MaterialPropertiesTable* aMaterialPropertiesTable =
	aMaterial->GetMaterialPropertiesTable();

      if (aMaterialPropertiesTable) {

	G4MaterialPropertyVector* theRefractionIndexVector = 
	  aMaterialPropertiesTable->GetProperty("RINDEX");

	if (theRefractionIndexVector) {
		
	  // Retrieve the first refraction index in vector
	  // of (photon momentum, refraction index) pairs 

	  G4double currentRI = (*theRefractionIndexVector)[0];

	  if (currentRI > 1.0) {

	    // Create first (photon momentum, Cerenkov Integral)
	    // pair  

	    G4double currentPM = theRefractionIndexVector->Energy(0);
	    G4double currentCAI = 0.0;

	    aPhysicsOrderedFreeVector->
	      InsertValues(currentPM , currentCAI);

	    // Set previous values to current ones prior to loop

	    G4double prevPM  = currentPM;
	    G4double prevCAI = currentCAI;
	    G4double prevRI  = currentRI;

	    // loop over all (photon momentum, refraction index)
	    // pairs stored for this material  

	    for (size_t i = 1; i < theRefractionIndexVector->GetVectorLength();i++)
	      {
		currentRI=(*theRefractionIndexVector)[i];
		currentPM = theRefractionIndexVector->Energy(i);

		currentCAI = 0.5*(1.0/(prevRI*prevRI) +
				  1.0/(currentRI*currentRI));

		currentCAI = prevCAI + 
		  (currentPM - prevPM) * currentCAI;

		aPhysicsOrderedFreeVector->
		  InsertValues(currentPM, currentCAI);

		prevPM  = currentPM;
		prevCAI = currentCAI;
		prevRI  = currentRI;
	      }

	  }
	}
      }

      // The Cerenkov integral for a given material
      // will be inserted in thePhysicsTable
      // according to the position of the material in
      // the material table. 

      thePhysicsTable->insertAt(i,aPhysicsOrderedFreeVector); 

    }
}

// GetMeanFreePath
// ---------------
//

G4double KM3Cherenkov::GetMeanFreePath(const G4Track&,
                                           G4double,
                                           G4ForceCondition*)
{
        return 1.;
}

G4double KM3Cherenkov::PostStepGetPhysicalInteractionLength(
                                           const G4Track& aTrack,
                                           G4double,
                                           G4ForceCondition* condition)
{
        *condition = NotForced;
        G4double StepLimit = DBL_MAX;

        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
        const G4Material* aMaterial = aTrack.GetMaterial();
        const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();

        G4double kineticEnergy = aParticle->GetKineticEnergy();
        const G4ParticleDefinition* particleType = aParticle->GetDefinition();
        G4double mass = particleType->GetPDGMass();

        // particle beta
        G4double beta = aParticle->GetTotalMomentum() /
	                aParticle->GetTotalEnergy();
        // particle gamma
        G4double gamma = aParticle->GetTotalEnergy()/mass;

        G4MaterialPropertiesTable* aMaterialPropertiesTable =
                            aMaterial->GetMaterialPropertiesTable();

        G4MaterialPropertyVector* Rindex = NULL;

        if (aMaterialPropertiesTable)
                     Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");

        G4double nMax;
        if (Rindex) {
           nMax = Rindex->GetMaxValue();
        } else {
           return StepLimit;
        }

        G4double BetaMin = 1./nMax;
        if ( BetaMin >= 1. ) return StepLimit;

        G4double GammaMin = 1./std::sqrt(1.-BetaMin*BetaMin);

        if (gamma < GammaMin ) return StepLimit;

        G4double kinEmin = mass*(GammaMin-1.);

        G4double RangeMin = G4LossTableManager::Instance()->
                                                   GetRange(particleType,
                                                            kinEmin,
                                                            couple);
        G4double Range    = G4LossTableManager::Instance()->
                                                   GetRange(particleType,
                                                            kineticEnergy,
                                                            couple);

        G4double Step = Range - RangeMin;
        if (Step < 1.*um ) return StepLimit;

        if (Step > 0. && Step < StepLimit) StepLimit = Step; 

        // If user has defined an average maximum number of photons to
        // be generated in a Step, then calculate the Step length for
        // that number of photons. 
 
        if (fMaxPhotons > 0) {

           // particle charge
           const G4double charge = aParticle->
                                   GetDefinition()->GetPDGCharge();

	   G4double MeanNumberOfPhotons = 
                    GetAverageNumberOfPhotons(charge,beta,aMaterial,Rindex);

           Step = 0.;
           if (MeanNumberOfPhotons > 0.0) Step = fMaxPhotons /
                                                 MeanNumberOfPhotons;

           if (Step > 0. && Step < StepLimit) StepLimit = Step;
        }

        // If user has defined an maximum allowed change in beta per step
        if (fMaxBetaChange > 0.) {

           G4double dedx = G4LossTableManager::Instance()->
                                                   GetDEDX(particleType,
                                                           kineticEnergy,
                                                           couple);

           G4double deltaGamma = gamma - 
                                 1./std::sqrt(1.-beta*beta*
                                                 (1.-fMaxBetaChange)*
                                                 (1.-fMaxBetaChange));

           Step = mass * deltaGamma / dedx;

           if (Step > 0. && Step < StepLimit) StepLimit = Step;

        }

        *condition = StronglyForced;
        return StepLimit;
}


// GetAverageNumberOfPhotons
// -------------------------
// This routine computes the number of Cerenkov photons produced per
// GEANT-unit (millimeter) in the current medium. 
//             ^^^^^^^^^^

G4double 
KM3Cherenkov::GetAverageNumberOfPhotons(const G4double charge,
					const G4double beta,
					const G4Material* aMaterial,
					G4MaterialPropertyVector* Rindex) const
{
  const G4double Rfact = 369.81/(eV * cm);

  if(beta <= 0.0)return 0.0;

  G4double BetaInverse = 1./beta;

  // Vectors used in computation of Cerenkov Angle Integral:
  // 	- Refraction Indices for the current material
  //	- new G4PhysicsOrderedFreeVector allocated to hold CAI's
 
  G4int materialIndex = aMaterial->GetIndex();

  // Retrieve the Cerenkov Angle Integrals for this material  

  G4PhysicsOrderedFreeVector* CerenkovAngleIntegrals =
    (G4PhysicsOrderedFreeVector*)((*thePhysicsTable)(materialIndex));

  if(!(CerenkovAngleIntegrals->IsFilledVectorExist()))return 0.0;

  // Min and Max photon momenta  
  G4double Pmin = Rindex->GetMinLowEdgeEnergy();
  G4double Pmax = Rindex->GetMaxLowEdgeEnergy();


  // Min and Max Refraction Indices 
  G4double nMin = Rindex->GetMinValue();	
  G4double nMax = Rindex->GetMaxValue();

  // Max Cerenkov Angle Integral 
  G4double CAImax = CerenkovAngleIntegrals->GetMaxValue();

  G4double dp, ge;

  // If n(Pmax) < 1/Beta -- no photons generated 

  if (nMax < BetaInverse) {
    //	  G4cout<<aParticle->GetTotalEnergy()<<" "<<aParticle->GetKineticEnergy()<<G4endl;
    //		dp = 0;
    //	ge = 0;
    return 0.0;
  } 

  // otherwise if n(Pmin) >= 1/Beta -- photons generated  

  else if (nMin > BetaInverse) {
    dp = Pmax - Pmin;	
    ge = CAImax; 
  } 

  // If n(Pmin) < 1/Beta, and n(Pmax) >= 1/Beta, then
  // we need to find a P such that the value of n(P) == 1/Beta.
  // Interpolation is performed by the GetPhotonMomentum() and
  // GetProperty() methods of the G4MaterialPropertiesTable and
  // the GetValue() method of G4PhysicsVector.  

  else {
    Pmin = Rindex->GetEnergy(BetaInverse);
    dp = Pmax - Pmin;

    // need boolean for current implementation of G4PhysicsVector
    // ==> being phased out
    G4bool isOutRange;
    G4double CAImin = CerenkovAngleIntegrals->
      GetValue(Pmin, isOutRange);
    ge = CAImax - CAImin;

    if (verboseLevel>0) {
      G4cout << "CAImin = " << CAImin << G4endl;
      G4cout << "ge = " << ge << G4endl;
    }
  }
	
  // Calculate number of photons 
  G4double NumPhotons = Rfact * charge/eplus * charge/eplus *
    (dp - ge * BetaInverse*BetaInverse);

  return NumPhotons;		
}

