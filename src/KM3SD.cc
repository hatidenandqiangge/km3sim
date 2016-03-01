//
#include "KM3SD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#ifdef G4TRACK_INFORMATION
#include "KM3TrackInformation.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KM3SD::KM3SD(std::string name) : G4VSensitiveDetector(name) {
  std::string HCname;
  collectionName.insert(HCname = "MyCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KM3SD::~KM3SD() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KM3SD::Initialize(G4HCofThisEvent *HCE) {
  MyCollection =
      new KM3HitsCollection(SensitiveDetectorName, collectionName[0]);
#ifdef G4MYFIT_PARAMETERIZATION
  int TotalNumberOfCathods = myStDetector->allCathods->GetNumberOfCathods();
  if (TotalNumberOfCathods > 20000)
    G4Exception("KM3SD::Initialize Number of cathods for energy fit is greater "
                "than 20000. Change the corresponding dimension in KM3SD.hh",
                "", FatalException, "");
  for (int i = 0; i < TotalNumberOfCathods; i++)
    ArrayParam[i] = 0;
  for (int i = 0; i < TotalNumberOfCathods; i++)
    ArrayParamAll[i] = 0;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool KM3SD::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist) {
  //  G4cout<< aStep->GetTrack()->GetDefinition()->GetParticleName()<<G4endl;
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() !=
      "opticalphoton") // this may have to change (do not kiil every particle on
                       // the photocathod)
  {
    //    aStep->GetTrack()->SetTrackStatus(fStopAndKill);   //kill every
    //    particle except photons that hit the cathod or dead part of OM
    //    G4cout<<"STOPPED"<<" particle "<<
    //    aStep->GetTrack()->GetDefinition()->GetParticleName() <<G4endl;
    return false;
  }

  // next kill photons that are incident on the back end of the pmt. The
  // definition of deadvolume is obsolete. so the following never happens
  //  if
  //  (aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName()=="DeadVolume_PV")
  //  //for newgeant  add "_PV" at the end of physical volume name
  //  {
  //    aStep->GetTrack()->SetTrackStatus(fStopAndKill); //kill a photon that
  //    hits the dead part of the OM
  //    //G4cout<<"STOPPED"<<G4endl;
  //    return false;
  //  }

  double edep = aStep->GetTrack()->GetTotalEnergy();
  if (edep == 0.) {
    aStep->GetTrack()->SetTrackStatus(
        fStopAndKill); // kill a photon with zero energy (?)
    return false;
  }

// comment about the fit,em,ha parametrizations
// In case that the Cathod is composed by water, there is a possibility that the
// photon
// is scatered inside the Cathod. In this case the hit is double in this Cathod.
// When it enters
// and when it leaves. In the setup with ~105000 cathods this increases the hit
// count by 0.2% when we have
// 60m scaterin length and 0.6% when we have 20m scattering length. This is
// negligible and I dont
// take this into account.
#ifdef G4MYFIT_PARAMETERIZATION
  /////next is new cathod id finding mode/////
  int Depth = aStep->GetPreStepPoint()->GetTouchable()->GetHistoryDepth();
  int History[10];
  for (int idep = 0; idep < Depth; idep++) {
    History[idep] = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(
        Depth - 1 - idep);
  }
  int id = myStDetector->allCathods->GetCathodId(Depth, History);
  G4ThreeVector posPMT = myStDetector->allCathods->GetPosition();
  double time = aStep->GetPostStepPoint()->GetGlobalTime();
  G4ThreeVector dirPARAM(0.0, 0.0, 1.0);
  double startz = -600.0 * meter;
  G4ThreeVector vertexPARAM(0.0, 0.0, startz);
  double TRes = TResidual(time, posPMT, vertexPARAM, dirPARAM);
  if (fabs(TRes) < 20.0 * ns)
    ArrayParam[id]++;
  if ((TRes > -20.0 * ns) && (TRes < 100.0 * ns))
    ArrayParamAll[id]++;
//  G4cout<<"ForFit "<<posPMT[0]/m<<" "<<posPMT[1]/m<<" "<<posPMT[2]/m<<"
//  "<<time<<" "<<TRes<<G4endl;
#endif
#if (defined(G4MYEM_PARAMETERIZATION) || defined(G4MYHA_PARAMETERIZATION)) &&  \
    !defined(G4MYK40_PARAMETERIZATION) // newha
  G4ThreeVector Dir = aStep->GetTrack()->GetMomentumDirection();
  // newmie don't count photons that are primary particles and have not been
  // scattered
  if (aStep->GetTrack()->GetParentID() == 0) {
    G4ThreeVector DirIni = aStep->GetTrack()->GetVertexMomentumDirection();
    double deviation = Dir.dot(DirIni);
    if (deviation > 0.99999999)
      return false;
  }
  // newmie
  int Depth = aStep->GetPreStepPoint()->GetTouchable()->GetHistoryDepth();
  int History[10];
  for (int idep = 0; idep < Depth; idep++) {
    History[idep] = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(
        Depth - 1 - idep);
  }
  int idcathod = myStDetector->allCathods->GetCathodId(Depth, History);
  //  double CathodRadius=myStDetector->allCathods->GetCathodRadius();
  G4ThreeVector PosHit = aStep->GetPreStepPoint()->GetPosition();
  //  double
  //  theDist1=(aStep->GetPreStepPoint()->GetPosition()-PosPMT).mag()/cm;
  //  //tempotest
  //  double
  //  theDist2=(aStep->GetPostStepPoint()->GetPosition()-PosPMT).mag()/cm;
  //  //tempotest
  //  G4cout <<"theMyDist "<<aStep->GetTrack()->GetTrackID()<<" "<<id<<"
  //  "<<theDist1<<" "<<theDist2<<G4endl;//tempotest
  G4ThreeVector FromGeneToOM = PosHit; // - myStDetector->MyGenerator->position;
                                       // //position is always 0,0,0
  double dist = FromGeneToOM.mag();
  //  if(fabs(dist-CathodRadius)>1.0*mm)G4Exception("KM3SD::ProcessHits in Param
  //  distances differ","",FatalException,"");
  double thetime = aStep->GetPreStepPoint()->GetGlobalTime();
  double cosangle1 =
      (myStDetector->MyGenerator->direction).dot(FromGeneToOM) / dist; // D.z
  double cosangle2 = (1.0 / dist) * FromGeneToOM.dot(Dir); // d.z
  double cosangle3;
  if (cosangle1 <= -1.0 || cosangle1 >= 1.0 || cosangle2 <= -1.0 ||
      cosangle2 >= 1.0) {
    cosangle3 = -1.0 + 2.0 * G4UniformRand();
  } else {
    cosangle3 =
        -(myStDetector->MyGenerator->direction.dot(Dir) -
          cosangle2 * cosangle1) /
        sqrt((1.0 - cosangle2 * cosangle2) * (1.0 - cosangle1 * cosangle1));
  }
  //  if(!(cosangle3>=-1.1 && cosangle3<=1.1))G4cout<<"NaNNaN "<<cosangle3<<"
  //  "<<cosangle2<<" "<<cosangle1<<G4endl;

  if (cosangle3 >= 1.0)
    cosangle3 = 0.0;
  else if (cosangle3 <= -1.0)
    cosangle3 = pi;
  else
    cosangle3 = acos(cosangle3);

  static int ooo = 0;
  if (ooo == 0) {
    ooo = 1;
    G4Material *aMaterial = G4Material::GetMaterial("Cathod");
    double MaxQE = -1;
    double PhEneAtMaxQE;
    G4MaterialPropertyVector *aPropertyVector =
        aMaterial->GetMaterialPropertiesTable()->GetProperty("Q_EFF");
    for (size_t i = 0; i < aPropertyVector->GetVectorLength(); i++) {
      double ThisQE = (*aPropertyVector)[i];
      double ThisPhEne = aPropertyVector->Energy(i);
      if (ThisQE > MaxQE) {
        MaxQE = ThisQE;
        PhEneAtMaxQE = ThisPhEne;
      }
    }
    aMaterial = G4Material::GetMaterial("Water");
    G4MaterialPropertyVector *GroupVel =
        aMaterial->GetMaterialPropertiesTable()->GetProperty("GROUPVEL");
    thespeedmaxQE = GroupVel->Value(PhEneAtMaxQE); // coresponds to the maximum
                                                   // qe each time. This is the
                                                   // right one
    VirtualAbsVector =
        aMaterial->GetMaterialPropertiesTable()->GetProperty("ABSLENGTH");
    TrueAbsVector =
        aMaterial->GetMaterialPropertiesTable()->GetProperty("ABSLENGTH_TRUE");
  }

  double theFastTime = dist / thespeedmaxQE;
  thetime -= theFastTime;

  // here we calculate the weight of each detected photon, because as simulation
  // abslength we have put a virtual one
  double tracklength = aStep->GetTrack()->GetTrackLength();
  double abs_virtual = VirtualAbsVector->Value(edep);
  double abs_true = TrueAbsVector->Value(edep);
  double weight = exp(-tracklength * (1.0 / abs_true - 1.0 / abs_virtual));
  //-----------------//

  //  if(dist>10*m)fprintf(myStDetector->outfile,"Param %12.4e %12.4e %12.4e
  //  %12.4e %12.4e %12.4e\n",  //tempo
  //  		       dist/m,cosangle1,cosangle2,cosangle3/degree,thetime,weight);
  //  //tempo
  // note: The binning in cosangle2,cosangle3 is not good enough
  // it must be done again
  // look folder ReBinc2c3 in Work/bin/Linux and below

  int distbin = idcathod; // the distances definition are recorded in the gdml
                            // geometry file
  // if(dist<6*meter){
  //   distbin=int(dist/meter);} //per 1 m up to 6m.
  // else if(dist<20*meter){ //per 2 meter from 6m to 20m
  //   distbin=int((dist-6.0*meter)/(meter*2.0))+6;}
  // else if(dist<50*meter){ //per 5 meter from 20 to 50m
  //   distbin=int((dist-20.0*meter)/(meter*5.0))+13;}
  // else if(dist<210*meter){ //per 20 meter from 50 to 210m
  //   distbin=int((dist-50.0*meter)/(meter*20.0))+19;}
  // else if(dist<510*meter){ //per 50 meter from 210 to 510m
  //   distbin=int((dist-210.0*meter)/(meter*50.0))+27;}
  // else distbin=32; //accumulate all above 510m

  // next we adjust the weight for the incident angle. Later in EventAction we
  // adjust it for normal surface
  // but first make sure that cosangle2 is not too small
  if (fabs(cosangle2) < 0.001) {
    if (cosangle2 > 0.0)
      cosangle2 = 0.001;
    else
      cosangle2 = -0.001;
  }
  weight /= fabs(cosangle2);
  //  if(fabs(cosangle2)<1.e-2)G4cout << "ppproblemmm "<<dist/m<<"
  //  "<<cosangle1<<" "<<cosangle2<<" "<<cosangle3/degree<<" "<<thetime<<G4endl;
  int cang1bin; // we estimate the theta1 bin
  // theta bin limits are:
  //-1 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 //15
  //bins
  // 0.500 0.525 0.550 0.575 0.600 0.625 0.650 0.675 0.700 //8 bins
  // 0.700 0.705 0.710 0.715 0.720 0.725 0.730 0.735 0.740 0.745 0.750 0.755
  // 0.760 0.765 0.770 0.775 0.780 0.785 0.790 0.795 0.800 //20 bins
  // 0.800 0.825 0.850 0.875 0.900 0.925 0.950 0.975 1.0 //8 bins

  // new theta binning (only for direct)
  //-1 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 //15
  //bins
  // 0.500 0.525 0.550 0.575 0.600 0.625 0.650 0.675 0.700 //8 bins
  // 0.700 0.705 0.710 0.715 0.720 0.725 //5 bins
  // 0.725 0.726 0.727 0.728 0.729 0.730 0.731 0.732 0.733 0.734 0.735 0.736
  // 0.737 0.738 0.739 0.740 0.741 0.742 0.743 0.744 0.745 0.746 0.747 0.748
  // 0.749 0.750 //25 bins
  // 0.750 0.755 0.760 0.765 0.770 0.775 0.780 0.785 0.790 0.795 0.800 //10 bins
  // 0.800 0.825 0.850 0.875 0.900 0.925 0.950 0.975 1.0 //8 bins

  bool FineBin = false;
  int VertexSolidAngleBins = 51;
  int TimeSolidAngleBins = 52;
  int TimeBins = 111;
  int OMSolidAngleBins = 834;
  if (myStDetector->MyGenerator->ParamEnergy == 0.0) {
    FineBin = true;
    VertexSolidAngleBins = 71;
  }

  if (!FineBin) {
    if (cosangle1 < 0.5) {
      cang1bin = int((cosangle1 + 1.0) / 0.1);
    } else if (cosangle1 < 0.7) {
      cang1bin = int((cosangle1 - 0.5) / 0.025) + 15;
    } else if (cosangle1 < 0.8) {
      cang1bin = int((cosangle1 - 0.7) / 0.005) + 23;
    } else if (cosangle1 < 1.0) {
      cang1bin = int((cosangle1 - 0.8) / 0.025) + 43;
    } else {
      cang1bin = VertexSolidAngleBins - 1;
    }
  } else {
    if (cosangle1 < 0.5) {
      cang1bin = int((cosangle1 + 1.0) / 0.1);
    } else if (cosangle1 < 0.7) {
      cang1bin = int((cosangle1 - 0.5) / 0.025) + 15;
    } else if (cosangle1 < 0.725) {
      cang1bin = int((cosangle1 - 0.7) / 0.005) + 23;
    } else if (cosangle1 < 0.75) {
      cang1bin = int((cosangle1 - 0.725) / 0.001) + 28;
    } else if (cosangle1 < 0.8) {
      cang1bin = int((cosangle1 - 0.75) / 0.005) + 53;
    } else if (cosangle1 < 1.0) {
      cang1bin = int((cosangle1 - 0.8) / 0.025) + 63;
    } else {
      cang1bin = VertexSolidAngleBins - 1;
    }
  }

  // change to be done
  // the variance of the distribution vs cosangle3 is due only to the
  // no zero length of the shower
  // this is zero for direct cherenkov, and almost zero for delta rays
  // however for em showers is not zero
  // better binning is done for the angle2 and angle3
  // 18 bins in cosangle2 and variable bins in cosangle3
  // the following:

  if (cosangle2 >= 1.0)
    cosangle2 = 0.0;
  else if (cosangle2 <= -1.0)
    cosangle2 = pi;
  else
    cosangle2 = acos(cosangle2);

  cosangle2 /= degree;
  cosangle3 /= degree;
  int cang23bin;
  if (cosangle2 < 4.43924) {
    if (cosangle3 < 55.0)
      cang23bin = 0;
    else if (cosangle3 < 105.0)
      cang23bin = 1;
    else
      cang23bin = 2;
  } else if (cosangle2 < 6.27957) {
    if (cosangle3 < 45.0)
      cang23bin = 3;
    else if (cosangle3 < 110.0)
      cang23bin = 4;
    else
      cang23bin = 5;
  } else if (cosangle2 < 8.10961) {
    if (cosangle3 < 45.0)
      cang23bin = 6;
    else if (cosangle3 < 110.0)
      cang23bin = 7;
    else
      cang23bin = 8;
  } else if (cosangle2 < 9.93636) {
    if (cosangle3 < 45.0)
      cang23bin = 9;
    else if (cosangle3 < 95.0)
      cang23bin = 10;
    else
      cang23bin = 11;
  } else if (cosangle2 < 11.4783) {
    if (cosangle3 < 45.0)
      cang23bin = 12;
    else if (cosangle3 < 120.0)
      cang23bin = 13;
    else
      cang23bin = 14;
  } else if (cosangle2 < 14.0699) {
    if (cosangle3 < 45.0)
      cang23bin = 15;
    else if (cosangle3 < 110.0)
      cang23bin = 16;
    else
      cang23bin = 17;
  } else if (cosangle2 < 18.1949) {
    if (cosangle3 < 40.0)
      cang23bin = 18;
    else if (cosangle3 < 75.0)
      cang23bin = 19;
    else
      cang23bin = 20;
  } else if (cosangle2 < 22.3316) {
    if (cosangle3 < 30.0)
      cang23bin = 21;
    else if (cosangle3 < 45.0)
      cang23bin = 22;
    else if (cosangle3 < 95.0)
      cang23bin = 23;
    else
      cang23bin = 24;
  } else if (cosangle2 < 25.8419) {
    if (cosangle3 < 40.0)
      cang23bin = 25;
    else if (cosangle3 < 75.0)
      cang23bin = 26;
    else if (cosangle3 < 110.0)
      cang23bin = 27;
    else
      cang23bin = 28;
  } else if (cosangle2 < 36.8699) {
    if (cosangle3 < 30.0)
      cang23bin = 29;
    else if (cosangle3 < 45.0)
      cang23bin = 30;
    else if (cosangle3 < 85.0)
      cang23bin = 31;
    else if (cosangle3 < 145.0)
      cang23bin = 32;
    else
      cang23bin = 33;
  } else if (cosangle2 < 45.5730) {
    if (cosangle3 < 40.0)
      cang23bin = 34;
    else if (cosangle3 < 75.0)
      cang23bin = 35;
    else if (cosangle3 < 115.0)
      cang23bin = 36;
    else
      cang23bin = 37;
  } else if (cosangle2 < 53.1301) {
    if (cosangle3 < 40.0)
      cang23bin = 38;
    else if (cosangle3 < 95.0)
      cang23bin = 39;
    else
      cang23bin = 40;
  } else if (cosangle2 < 63.2563) {
    if (cosangle3 < 45.0)
      cang23bin = 41;
    else if (cosangle3 < 85.0)
      cang23bin = 42;
    else
      cang23bin = 43;
  } else if (cosangle2 < 90.0000) {
    if (cosangle3 < 40.0)
      cang23bin = 44;
    else if (cosangle3 < 110.0)
      cang23bin = 45;
    else
      cang23bin = 46;
  } else if (cosangle2 < 113.578) {
    if (cosangle3 < 65.0)
      cang23bin = 47;
    else if (cosangle3 < 120.0)
      cang23bin = 48;
    else
      cang23bin = 49;
  } else if (cosangle2 < 143.130) {
    cang23bin = 50;
  } else {
    cang23bin = TimeSolidAngleBins - 1;
  }

  static double dtheta = -1;
  if (dtheta < 0) {
    // definition of two solid angle areas, one with 3 degrees binning and the
    // second with 6 degrees binning
    int NumberOfThetas1 = 34;
    int NumberOfThetas2 = 13;
    double Theta1Min = 0.0;
    double Theta1Max = 102.0;
    double Theta2Min = 102.0;
    double Theta2Max = 180.0;

    NumberOfThetas = NumberOfThetas1 + NumberOfThetas2;
    double ibinNum_Tot = 0;
    // first solid angle area
    dtheta = (Theta1Max - Theta1Min) / NumberOfThetas1;
    double cosdtheta = cos(dtheta * degree);
    for (int ith = 0; ith < NumberOfThetas1; ith++) {
      double thetalow = Theta1Min + dtheta * ith;
      double thetahigh = thetalow + dtheta;
      theta_Low[ith] = thetalow;
      theta_High[ith] = thetahigh;
      double costhetalow = fabs(cos(thetalow * degree));
      double costhetahigh = fabs(cos(thetahigh * degree));
      double cosmin;
      if (costhetalow < costhetahigh)
        cosmin = costhetalow;
      else
        cosmin = costhetahigh;
      cosmin = cosmin * cosmin;
      double cosdphi = (cosdtheta - cosmin) / (1 - cosmin);
      double dphi = acos(cosdphi) / degree;
      if (dphi < 9.0)
        dphi = 9.0;
      NumberOfPhis[ith] = int(ceil(180.0 / dphi));
      dphi = 180.0 / NumberOfPhis[ith];
      for (int iph = 0; iph < NumberOfPhis[ith]; iph++) {
        double philow = dphi * iph;
        double phihigh = philow + dphi;
        phi_Low[ith][iph] = philow;
        phi_High[ith][iph] = phihigh;
        ibinNum_Tot++;
      }
      phi_High[ith][NumberOfPhis[ith] - 1] = 180.0;
    }
    theta_High[NumberOfThetas1 - 1] = Theta1Max;
    // second sold angle area
    dtheta = (Theta2Max - Theta2Min) / NumberOfThetas2;
    cosdtheta = cos(dtheta * degree);
    for (int ith = NumberOfThetas1; ith < NumberOfThetas; ith++) {
      double thetalow = Theta2Min + dtheta * (ith - NumberOfThetas1);
      double thetahigh = thetalow + dtheta;
      theta_Low[ith] = thetalow;
      theta_High[ith] = thetahigh;
      double costhetalow = fabs(cos(thetalow * degree));
      double costhetahigh = fabs(cos(thetahigh * degree));
      double cosmin;
      if (costhetalow < costhetahigh)
        cosmin = costhetalow;
      else
        cosmin = costhetahigh;
      cosmin = cosmin * cosmin;
      double cosdphi = (cosdtheta - cosmin) / (1 - cosmin);
      double dphi = acos(cosdphi) / degree;
      if (dphi < 9.0)
        dphi = 9.0;
      NumberOfPhis[ith] = int(ceil(180.0 / dphi));
      dphi = 180.0 / NumberOfPhis[ith];
      for (int iph = 0; iph < NumberOfPhis[ith]; iph++) {
        double philow = dphi * iph;
        double phihigh = philow + dphi;
        phi_Low[ith][iph] = philow;
        phi_High[ith][iph] = phihigh;
        ibinNum_Tot++;
      }
      phi_High[ith][NumberOfPhis[ith] - 1] = 180.0;
    }
    theta_High[NumberOfThetas - 1] = Theta2Max;
    /////

    G4cout << "------Summary of solid angle splitting for calculation of "
              "photon numbers------"
           << G4endl;
    G4cout << "NumberOfThetas= " << NumberOfThetas
           << " Total Number of bins= " << ibinNum_Tot << G4endl;
    int count = 0;
    for (int ith = 0; ith < NumberOfThetas; ith++) {
      G4cout << "      NumberOfPhis= " << NumberOfPhis[ith]
             << " ThetaLimits= " << theta_Low[ith] << " - " << theta_High[ith]
             << G4endl;
      for (int iph = 0; iph < NumberOfPhis[ith]; iph++) {
        G4cout << count << "          PhiLimits= " << phi_Low[ith][iph] << " - "
               << phi_High[ith][iph] << G4endl;
        count++;
      }
    }
    G4cout << "------Summary of solid angle splitting for calculation of "
              "photon numbers------"
           << G4endl;
    if (ibinNum_Tot != OMSolidAngleBins)
      G4Exception("KM3SD::ProcessHits Number of calculated bins is not equal "
                  "to OMSolidAngleBins",
                  "", FatalException, "");
  }

  int ibinNum23 = 0;
  int ith;
  for (ith = 0; ith < NumberOfThetas - 1; ith++) {
    if (cosangle2 >= theta_High[ith]) {
      for (int iph = 0; iph < NumberOfPhis[ith]; iph++)
        ibinNum23++;
    } else
      break;
  }
  for (int iph = 0; iph < NumberOfPhis[ith]; iph++) {
    if (cosangle3 >= phi_High[ith][iph]) {
      ibinNum23++;
    } else
      break;
  }
  if (cosangle3 >= 180.0)
    ibinNum23--;

  int ibint;
  if (thetime < 10) {
    ibint = int((thetime + 10.0) / 0.5);
    if (ibint < 0)
      ibint = 0;
  } else if (thetime < 30) {
    ibint = int(thetime - 10.0) + 40;
  } else if (thetime < 90) {
    ibint = int((thetime - 30.0) / 3.0) + 60;
  } else if (thetime < 200) {
    ibint = int((thetime - 90.0) / 11.0) + 80;
  } else if (thetime < 1000) {
    ibint = int((thetime - 200.0) / 50.0) + 90;
  } else if (thetime < 2000) {
    ibint = int((thetime - 1000.0) / 200.0) + 106;
  } else
    ibint = TimeBins - 1;

  int ibin_d1 = cang1bin + VertexSolidAngleBins * distbin;
  (*myPhotonsNumber)[ibin_d1] += (long double)weight;
  int ibin_d123 = cang23bin + TimeSolidAngleBins * ibin_d1;
  int ibin_d123t = ibint + TimeBins * ibin_d123;
  (*myPhotonsTime)[ibin_d123t] += (long double)weight;
  int ibin_d123num =
      ibinNum23 +
      OMSolidAngleBins * ibin_d1; // OMSolidAngleBins is ibinNum_Tot above
  (*myPhotonsTh2Th3Num)[ibin_d123num] += (long double)weight;
  (*myPhotonsTh2)[ibin_d123num] += (long double)(weight * cosangle2);
  (*myPhotonsTh3)[ibin_d123num] += (long double)(weight * cosangle3);

#endif
#ifndef G4MYFIT_PARAMETERIZATION
#if (!defined(G4MYEM_PARAMETERIZATION) &&                                      \
     !defined(G4MYHA_PARAMETERIZATION)) ||                                     \
    defined(G4MYK40_PARAMETERIZATION) // newha
  if (MyCollection->entries() < 10000000) {
    G4ThreeVector photonDirection = aStep->GetTrack()->GetMomentumDirection();

// newmie
// here we disgard photons that have been scattered and are created with
// parametrization also in KM3Cherenkov
// this for parametrization running. Not anymore, since these are killed in the
// G4OpMie scattering process
#ifdef G4TRACK_INFORMATION
    KM3TrackInformation *info = NULL;
#endif
    //#ifndef G4DISABLE_PARAMETRIZATION
    //#ifdef G4TRACK_INFORMATION
    //    info= (KM3TrackInformation*)(aStep->GetTrack()->GetUserInformation());
    //    if(info->GetEmittedAsScattered()){
    //      G4ThreeVector DirIni=
    //      aStep->GetTrack()->GetVertexMomentumDirection();
    //      double deviation=photonDirection.dot(DirIni);
    //      if(deviation<0.99999999){
    //	aStep->GetTrack()->SetTrackStatus(fStopAndKill); //it is killed, because
    //it is useless
    //	return false;
    //      }
    //    }
    //#else
    //    G4Exception("KM3SD::ProcessHits: Parametrization application needs
    //    track information","",FatalException,"");
    //#endif
    //#endif
    // newmie

    /////next is new cathod id finding mode/////
    int Depth = aStep->GetPreStepPoint()->GetTouchable()->GetHistoryDepth();
    int History[10];
    for (int idep = 0; idep < Depth; idep++) {
      History[idep] =
          aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(Depth - 1 -
                                                                     idep);
    }
    int id = myStDetector->allCathods->GetCathodId(Depth, History);
    //////check if this photon passes after the angular acceptance//////
    G4ThreeVector PMTDirection = myStDetector->allCathods->GetDirection();
    double CathodRadius = myStDetector->allCathods->GetCathodRadius();
    double CathodHeight = myStDetector->allCathods->GetCathodHeight();
    if (!AcceptAngle(photonDirection.dot(PMTDirection), CathodRadius,
                     CathodHeight, false)) {
      // at this point we dont kill the track if it is not accepted due to
      // anglular acceptance
      // this has an observable effect a few percent only when running
      // simulation with
      // parametrization turned off and sparce detectors.
      // However with dense detectors ~0.04 OMs/m^3 there is a 30% effect, and
      // 0.005OMs/m^3 a 3.5% effect, based on the covered solid angle
      return false;
    }
    //    G4cout<<"Accepted"<<G4endl;
    //    G4cout<<"----------------------------------------------------------------"<<G4endl;

    KM3Hit *newHit = new KM3Hit();
    newHit->SetCathodId(id);
    newHit->SetTime(aStep->GetPostStepPoint()->GetGlobalTime());

    int originalTrackCreatorProcess;
    int originalParentID;
#ifdef G4TRACK_INFORMATION
    if (info == NULL)
      info = (KM3TrackInformation *)(aStep->GetTrack()->GetUserInformation());
    originalParentID = info->GetOriginalParentID();
    std::string creator = info->GetOriginalTrackCreatorProcess();
    if (creator == "KM3Cherenkov")
      originalTrackCreatorProcess = 0;
    else if (creator == "muPairProd")
      originalTrackCreatorProcess = 1;
    else if (creator == "muIoni")
      originalTrackCreatorProcess = 2;
    else if (creator == "muBrems")
      originalTrackCreatorProcess = 3;
    else if (creator == "muonNuclear")
      originalTrackCreatorProcess = 4;
    else if (creator == "Decay")
      originalTrackCreatorProcess = 8;
    else if (creator == "muMinusCaptureAtRest")
      originalTrackCreatorProcess = 9;
    else
      originalTrackCreatorProcess = 5;
#ifdef G4MYLASER_PARAMETERIZATION
    info->KeepScatteringPosition(aStep->GetPostStepPoint()->GetPosition(), 1.0);
    int NumberOfScatters = info->GetNumberOfScatters();
    newHit->SetIManyScatters(NumberOfScatters - 1);
    for (int is = 0; is < NumberOfScatters - 1; is++) {
      newHit->SetScatteringSteps(is, (info->GetScatteringPosition(is + 1) -
                                      info->GetScatteringPosition(is))
                                         .mag());
      newHit->SetScatteringAngles(is, info->GetScatteringAngle(is + 1));
    }
#endif
#else
    originalParentID = 1;
    originalTrackCreatorProcess = 0;
#endif
    int originalInfo;
    originalInfo = (originalParentID - 1) * 10 + originalTrackCreatorProcess;
    //    newHit->SetoriginalInfo(int(1.e6*h_Planck*c_light/aStep->GetTrack()->GetTotalEnergy()));
    newHit->SetoriginalInfo(originalInfo);
    newHit->SetMany(1);

    // short    G4ThreeVector posHit=aStep->GetPostStepPoint()->GetPosition();
    // short    G4ThreeVector posPMT=myStDetector->allCathods->GetPosition();
    // short    G4ThreeVector posRel=posHit-posPMT;

    // short    double angleThetaIncident,anglePhiIncident;
    // short    double angleThetaDirection,anglePhiDirection;

    // short    angleThetaIncident=posRel.theta();
    // short    anglePhiIncident=posRel.phi();
    // short    angleThetaDirection=photonDirection.theta();
    // short    anglePhiDirection=photonDirection.phi();
    // convert to degrees
    // short    angleThetaIncident *= 180./M_PI;
    // short    anglePhiIncident *= 180./M_PI;
    // short    angleThetaDirection *= 180./M_PI;
    // short    anglePhiDirection *= 180./M_PI;
    // short    if(anglePhiIncident < 0.0)anglePhiIncident += 360.0;
    // short    if(anglePhiDirection < 0.0)anglePhiDirection += 360.0;

    // short    int angleIncident,angleDirection;
    // short    angleIncident = (int)(nearbyint(angleThetaIncident)*1000.0 +
    // nearbyint(anglePhiIncident));
    // short    angleDirection = (int)(nearbyint(angleThetaDirection)*1000.0 +
    // nearbyint(anglePhiDirection));

    // short    newHit->SetangleIncident(angleIncident);
    // short    newHit->SetangleDirection(angleDirection);

    MyCollection->insert(newHit);
  }
#endif
#endif

#if (!defined(G4MYEM_PARAMETERIZATION) &&                                      \
     !defined(G4MYHA_PARAMETERIZATION)) ||                                     \
    defined(G4MYK40_PARAMETERIZATION) // newha
#ifndef G4MYFIT_PARAMETERIZATION
  // killing must not been done, when we have EM or HA or FIT parametrizations
  // but it must be done for normal run, especially K40, SN and laser since
  // coincidences play big role there
  aStep->GetTrack()->SetTrackStatus(fStopAndKill); // kill the detected photon
#endif
#endif

  return true;
}
// this method is used to add hits from the EM shower model
// short void KM3SD::InsertExternalHit(int id,double time,int
// originalInfo,int angleDirection,int angleIncident)
void KM3SD::InsertExternalHit(int id, const G4ThreeVector &OMPosition,
                              double time, int originalInfo,
                              const G4ThreeVector &photonDirection) {
  /////////calculate the photon speed at max QE to correct time
  static int ooo = 0;
  if (ooo == 0) {
    ooo = 1;
    G4Material *aMaterial = G4Material::GetMaterial("Cathod");
    double MaxQE = -1;
    double PhEneAtMaxQE;
    G4MaterialPropertyVector *aPropertyVector =
        aMaterial->GetMaterialPropertiesTable()->GetProperty("Q_EFF");
    for (size_t i = 0; i < aPropertyVector->GetVectorLength(); i++) {
      double ThisQE = (*aPropertyVector)[i];
      double ThisPhEne = aPropertyVector->Energy(i);
      if (ThisQE > MaxQE) {
        MaxQE = ThisQE;
        PhEneAtMaxQE = ThisPhEne;
      }
    }
    aMaterial = G4Material::GetMaterial("Water");
    G4MaterialPropertyVector *GroupVel =
        aMaterial->GetMaterialPropertiesTable()->GetProperty("GROUPVEL");
    thespeedmaxQE = GroupVel->Value(PhEneAtMaxQE); // coresponds to the maximum
                                                   // qe each time. This is the
                                                   // right one
  }
//////////////////////////////////////////////////////////////
#ifndef G4MYFIT_PARAMETERIZATION
  if (MyCollection->entries() < 10000000) {
    G4ThreeVector PMTDirection = myStDetector->allCathods->GetDirection(id);
    //    G4cout << "OutFromParam "<<id<<"
    //    "<<photonDirection.dot(PMTDirection)<<G4endl;
    if (!AcceptAngle(photonDirection.dot(PMTDirection), 1.0, 1.0,
                     true)) { // the two 1.0 are the cathod radius and height
                              // that do not play any role in parametrization
      return; // it is not accepted
    }
    // correct the time to correspond to the cathod positions and not the OM
    // position
    G4ThreeVector PMTPosition = myStDetector->allCathods->GetPosition(id);
    double Tcorr =
        (photonDirection.dot(PMTPosition - OMPosition)) / thespeedmaxQE;
    KM3Hit *newHit = new KM3Hit();
    newHit->SetCathodId(id);
    newHit->SetTime(time + Tcorr);
    // short    newHit->SetangleIncident(angleIncident);
    // short    newHit->SetangleDirection(angleDirection);
    newHit->SetoriginalInfo(originalInfo);
    newHit->SetMany(1);
    MyCollection->insert(newHit);
  }
#else
  G4ThreeVector dirPARAM(0.0, 0.0, 1.0);
  double startz = -600.0 * meter;
  G4ThreeVector vertexPARAM(0.0, 0.0, startz);
  G4ThreeVector posPMT = myStDetector->allCathods->GetPosition(id);
  double TRes = TResidual(time, posPMT, vertexPARAM, dirPARAM);
  if (fabs(TRes) < 20.0 * ns)
    ArrayParam[id]++;
  if ((TRes > -20.0 * ns) && (TRes < 100.0 * ns))
    ArrayParamAll[id]++;
//  G4cout<<"ForFit "<<posPMT[0]/m<<" "<<posPMT[1]/m<<" "<<posPMT[2]/m<<"
//  "<<time<<" "<<TRes<<G4endl;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#if (!defined(G4MYEM_PARAMETERIZATION) &&                                      \
     !defined(G4MYHA_PARAMETERIZATION)) ||                                     \
    defined(G4MYK40_PARAMETERIZATION) // newha
void KM3SD::EndOfEvent(G4HCofThisEvent *HCE) {
#ifdef G4MYK40_PARAMETERIZATION
  // here get access to the primary generation action generation parameters
  // idbeam,ParamEnergy,random_R
  int idbeam = myStDetector->MyGenerator->idbeam;
  double ParamEnergy = myStDetector->MyGenerator->ParamEnergy / keV;
  double random_R = myStDetector->MyGenerator->random_R / meter;
#ifdef G4MYLASER_PARAMETERIZATION
  idbeam = int(idbeam / myStDetector->Quantum_Efficiency);
#endif
#endif

  if (verboseLevel > 0) {
    int TotalNumberOfCathods = myStDetector->allCathods->GetNumberOfCathods();
#ifdef G4MYFIT_PARAMETERIZATION
    outfile = myStDetector->outfilePar;
    double ParamMomentum =
        myStDetector->event_action->centerMomentum[0] / GeV;
    for (int ica = 0; ica < TotalNumberOfCathods; ica++) {
      G4ThreeVector PosPMT = myStDetector->allCathods->GetPosition(ica);
      double ParamDistance =
          sqrt(PosPMT[0] * PosPMT[0] + PosPMT[1] * PosPMT[1]) / meter;
      fprintf(outfile, "%.4e %.4e %d %d\n", ParamMomentum, ParamDistance,
              ArrayParam[ica], ArrayParamAll[ica]);
    }
#else
    outfile = myStDetector->outfile;
    int NbHits = MyCollection->entries();
    static int ooo = 0; // count total
    ooo += NbHits; // count total
    G4cout << "Total Hits: " << ooo << G4endl; // count total
    G4cout << "This Event Hits: " << NbHits << G4endl; // count for this event
    int i;

    // here we sort MyCollection according to ascending pmt number
    std::vector<KM3Hit *> *theCollectionVector = MyCollection->GetVector();
    QuickSort(0, theCollectionVector, 0, NbHits - 1);
    //-------//from now on the hits are sorted in cathod id
    // next sort according to ascending time for each cathod id
    if (NbHits > 1) {
      int prevCathod, currentCathod;
      int istart, istop;
      istart = 0;
      prevCathod = (*MyCollection)[istart]->GetCathodId();
      for (i = 1; i < NbHits; i++) {
        currentCathod = (*MyCollection)[i]->GetCathodId();
        if (currentCathod != prevCathod) {
          istop = i - 1;
          QuickSort(1, theCollectionVector, istart, istop);
#ifndef G4MYK40_PARAMETERIZATION
          double MergeWindow = 0.5 * ns;
          MergeHits(istart, istop + 1, MergeWindow);
#endif
#ifdef G4MYLASER_PARAMETERIZATION
          double MergeWindow = 0.1 * ns;
          MergeHits(istart, istop + 1, MergeWindow);
#endif
          prevCathod = currentCathod;
          istart = i;
        } else if ((currentCathod == prevCathod) && (i == NbHits - 1)) {
          istop = i;
          QuickSort(1, theCollectionVector, istart, istop);
#ifndef G4MYK40_PARAMETERIZATION
          double MergeWindow = 0.5 * ns;
          MergeHits(istart, istop + 1, MergeWindow);
#endif
#ifdef G4MYLASER_PARAMETERIZATION
          double MergeWindow = 0.1 * ns;
          MergeHits(istart, istop + 1, MergeWindow);
#endif
        }
      }
    }

    // find the number of hit entries to write
    int NbHitsWrite = 0;
    for (i = 0; i < NbHits; i++)
      if ((*MyCollection)[i]->GetMany() > 0)
        NbHitsWrite++;

    // find earliest hit time
    double timefirst = 1E20;
    for (i = 0; i < NbHits; i++) {
      if ((*MyCollection)[i]->GetTime() < timefirst &&
          (*MyCollection)[i]->GetMany() > 0)
        timefirst = (*MyCollection)[i]->GetTime();
    }

    // find how many cathods are hitted
    int allhit = 0;
    int prevcathod = -1;
    for (i = 0; i < NbHits; i++) {
      if (prevcathod != (*MyCollection)[i]->GetCathodId())
        allhit++;
      prevcathod = (*MyCollection)[i]->GetCathodId();
    }
#ifdef G4MYK40_PARAMETERIZATION
    if (NbHits > 0) {
#ifdef G4MYLASER_PARAMETERIZATION
      fprintf(outfile, "%d\n", idbeam);
#else
      fprintf(outfile, "%d %.6e %.6e\n", idbeam, ParamEnergy, random_R);
#endif
      fprintf(outfile, "%.7e %d\n", timefirst * 1E-9, allhit);
    }
#else
    if (!myStDetector->useANTARESformat)
      fprintf(outfile, "%.7e %d\n", timefirst * 1E-9, allhit);
    else
      myStDetector->TheEVTtoWrite->AddNumberOfHits(NbHitsWrite);
#endif

    // find the last pmt how many hits has
    int LastPmtNumber;
    int LastHitNumber;
    for (i = NbHits - 1; i >= 0; i--)
      if ((*MyCollection)[i]->GetMany() > 0) {
        LastPmtNumber = (*MyCollection)[i]->GetCathodId();
        LastHitNumber = i;
        break;
      }
    int LastPmtNumHits = 0;
    for (i = NbHits - 1; i >= 0; i--)
      if ((*MyCollection)[i]->GetMany() > 0) {
        if ((*MyCollection)[i]->GetCathodId() == LastPmtNumber)
          LastPmtNumHits++;
        else
          break;
      }

    // new
    int numphotons = 0;
    double firstphoton = 1.E50;
    int numpes = 0;
    if (NbHits > 0)
      prevcathod = (*MyCollection)[0]->GetCathodId();
    int prevstart = 0;
    int numhit = 0;
    for (i = 0; i < NbHits; i++) {
      if ((*MyCollection)[i]->GetMany() > 0) {
        if (prevcathod == (*MyCollection)[i]->GetCathodId()) {
          numphotons++;
          numpes += (*MyCollection)[i]->GetMany();
          if ((*MyCollection)[i]->GetTime() - timefirst < firstphoton)
            firstphoton = (*MyCollection)[i]->GetTime() - timefirst;
        } else {
          if (myStDetector->vrmlhits) { // draw hits
            G4ThreeVector Cposition =
                myStDetector->allCathods->GetPosition(prevcathod);
            DrawCathodHit(numpes, Cposition);
          }
          if (!myStDetector->useANTARESformat)
            fprintf(outfile, "%d %d %.7e\n", prevcathod, numphotons,
                    firstphoton * 1E-9);
          for (int j = prevstart; j < i; j++) {
            if ((*MyCollection)[j]->GetMany() > 0) {
              numhit++;
              if (!myStDetector->useANTARESformat) {
                fprintf(outfile, "%.7e %d %d\n",
                        ((*MyCollection)[j]->GetTime() - timefirst) * 1E-9,
                        (*MyCollection)[j]->GetoriginalInfo(),
                        (*MyCollection)[j]->GetMany());
#if defined(G4MYLASER_PARAMETERIZATION) && defined(G4TRACK_INFORMATION)
                int IManyScatters = (*MyCollection)[j]->GetIManyScatters();
                fprintf(outfile, "%d\n", IManyScatters);
                for (int is = 0; is < IManyScatters; is++)
                  fprintf(outfile, "%.3e %.6e\n",
                          (*MyCollection)[j]->GetScatteringSteps(is) / m,
                          (*MyCollection)[j]->GetScatteringAngles(is));
#endif
              } else {
                // here write antares format info
                int originalInfo = (*MyCollection)[j]->GetoriginalInfo();
                int originalParticleNumber = originalInfo / 10 + 1;
                int originalTrackCreatorProcess =
                    originalInfo - (originalParticleNumber - 1) * 10;
                myStDetector->TheEVTtoWrite->AddHit(
                    numhit, prevcathod, double((*MyCollection)[j]->GetMany()),
                    (*MyCollection)[j]->GetTime(), originalParticleNumber,
                    (*MyCollection)[j]->GetMany(),
                    (*MyCollection)[j]->GetTime(), originalTrackCreatorProcess);
              }
            }
          }
          prevstart = i;
          numphotons = 1;
          firstphoton = (*MyCollection)[i]->GetTime() - timefirst;
          numpes = (*MyCollection)[i]->GetMany();
        }
        prevcathod = (*MyCollection)[i]->GetCathodId();
        //	if(i == (NbHits-1) ){
        //	if(numhit == (NbHitsWrite-1) ){
        if (numhit == (NbHitsWrite - LastPmtNumHits) && i == LastHitNumber) {
          if (myStDetector->vrmlhits) { // draw hits
            G4ThreeVector Cposition = myStDetector->allCathods->GetPosition(
                (*MyCollection)[i]->GetCathodId());
            DrawCathodHit(numpes, Cposition);
          }
          if (!myStDetector->useANTARESformat)
            fprintf(outfile, "%d %d %.7e\n", (*MyCollection)[i]->GetCathodId(),
                    numphotons, firstphoton * 1E-9);
          for (int j = prevstart; j < NbHits; j++) {
            if ((*MyCollection)[j]->GetMany() > 0) {
              numhit++;
              if (!myStDetector->useANTARESformat) {
                fprintf(outfile, "%.7e %d %d\n",
                        ((*MyCollection)[j]->GetTime() - timefirst) * 1E-9,
                        (*MyCollection)[j]->GetoriginalInfo(),
                        (*MyCollection)[j]->GetMany());
#if defined(G4MYLASER_PARAMETERIZATION) && defined(G4TRACK_INFORMATION)
                int IManyScatters = (*MyCollection)[j]->GetIManyScatters();
                fprintf(outfile, "%d\n", IManyScatters);
                for (int is = 0; is < IManyScatters; is++)
                  fprintf(outfile, "%.3e %.6e\n",
                          (*MyCollection)[j]->GetScatteringSteps(is) / m,
                          (*MyCollection)[j]->GetScatteringAngles(is));
#endif
              } else {
                // here write antares format info
                int originalInfo = (*MyCollection)[j]->GetoriginalInfo();
                int originalParticleNumber = originalInfo / 10 + 1;
                int originalTrackCreatorProcess =
                    originalInfo - (originalParticleNumber - 1) * 10;
                myStDetector->TheEVTtoWrite->AddHit(
                    numhit, (*MyCollection)[i]->GetCathodId(),
                    double((*MyCollection)[j]->GetMany()),
                    (*MyCollection)[j]->GetTime(), originalParticleNumber,
                    (*MyCollection)[j]->GetMany(),
                    (*MyCollection)[j]->GetTime(), originalTrackCreatorProcess);
              }
            }
          }
        }
      }
    }

// old
// int numphotons;double firstphoton;
// for(int ica=0;ica<TotalNumberOfCathods;ica++){ //for all benthos
//   numphotons=0;firstphoton=1.E50;
//   for(i=0;i<NbHits;i++){
// 	if((*MyCollection)[i]->GetCathodId()==ica){
// 	  numphotons++;  // hit number for this benthos
// 	  if((*MyCollection)[i]->GetTime()-timefirst<firstphoton)firstphoton=(*MyCollection)[i]->GetTime()-timefirst;
// 	}
//   }
//   if(numphotons>0){
// 	if(myStDetector->vrmlhits){  //draw hits
// 	  G4ThreeVector Cposition=myStDetector->allCathods->GetPosition(ica);
// 	  DrawCathodHit(numphotons,Cposition);
// 	}
// 	fprintf(outfile,"%d %d %.7e
// %d\n",ica,numphotons,firstphoton*1E-9,numphotons);
// 	for (i=0;i<NbHits;i++){
// 	  if((*MyCollection)[i]->GetCathodId()==ica)
// 	    fprintf(outfile,"%.7e %d %d %d\n",
// 		    ((*MyCollection)[i]->GetTime()-timefirst)*1E-9,
// 		    (*MyCollection)[i]->GetangleIncident(),
// 		    (*MyCollection)[i]->GetangleDirection(),
// 		    (*MyCollection)[i]->GetoriginalInfo());
// 	}
//   }
// }
// old
#endif
    if (myStDetector->vrmlhits) {
      static int HCID = -1;
      if (HCID < 0) {
        HCID = GetCollectionID(0);
      }
      HCE->AddHitsCollection(HCID, MyCollection);
    } else
      delete MyCollection;

    // G4cout << "\n-------->Hits Collection: in this event they are " <<
    // MyCollection->entries() << G4endl;
    // for (int i=0;i<NbHits;i++) (*MyCollection)[i]->Print();
  }
}
#else
void KM3SD::EndOfEvent(G4HCofThisEvent *HCE) { delete MyCollection; }
#endif

void KM3SD::MergeHits(int nfirst, int nlast, double MergeWindow) {
  if (nlast - nfirst < 2)
    return;
  int iuu, imany;
  int istart = nfirst;
  int istop = nfirst;
go77:
  istart = istop + 1;
  for (iuu = istart + 1; iuu <= nlast; iuu++) {
    if (((*MyCollection)[iuu - 1]->GetTime() -
         (*MyCollection)[istart - 1]->GetTime()) > MergeWindow) {
      istop = iuu - 1;
      goto go78;
    }
  }
  istop = iuu - 1;
go78:
  if (istart > nlast)
    return;
  imany = istop - istart + 1;
  if (imany > 1) {
    double MeanTime = 0.0;
    for (iuu = istart; iuu <= istop; iuu++)
      MeanTime += (*MyCollection)[iuu - 1]->GetTime();
    MeanTime /= imany;
    (*MyCollection)[istart - 1]->SetTime(MeanTime);
    (*MyCollection)[istart - 1]->SetMany(imany);
    for (iuu = istart + 1; iuu <= istop; iuu++)
      (*MyCollection)[iuu - 1]->SetMany(0);
  }
  goto go77;
}

#include "G4VisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

// Hit Draw Method (colours are descibing the number of photons, blue->red)
void KM3SD::DrawCathodHit(int NumberOfPhotons, G4ThreeVector pos) {
  G4VVisManager *pVVisManager = G4VVisManager::GetConcreteInstance();
  if (pVVisManager) {
    G4Circle circle(pos);
    circle.SetWorldRadius(220.0);
    circle.SetFillStyle(G4Circle::filled);
    double nphotons = double(NumberOfPhotons);
    if (nphotons > 100.0)
      nphotons = 100.0;
    double redcol = log10(double(nphotons)) / 2.0;
    G4Colour colour(redcol, 0., 1.0 - redcol);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

int KM3SD::ProcessMyCollection(KM3HitsCollection *aCollection) { return (0); }

void KM3SD::clear() {}

void KM3SD::PrintAll() {}

double KM3SD::TResidual(double time, const G4ThreeVector &position,
                          const G4ThreeVector &vertex,
                          const G4ThreeVector &dir) {
  double tnthc = 0.961; // this value depends on qe and water properties
  double ag, bg;
  G4ThreeVector Hit = position - vertex;
  ag = dir.dot(Hit);
  Hit -= ag * dir;
  bg = Hit.mag();
  return time - (ag + bg * tnthc) / c_light;
}

// Quick Sort Functions for Ascending Order
void KM3SD::QuickSort(int shorttype,
                      std::vector<KM3Hit *> *theCollectionVector, int top,
                      int bottom) {
  // top = subscript of beginning of array
  // bottom = subscript of end of array

  int middle;
  if (top < bottom) {
    if (shorttype == 0)
      middle = partition_CathodId(theCollectionVector, top, bottom);
    else
      middle = partition_Time(theCollectionVector, top, bottom);
    QuickSort(shorttype, theCollectionVector, top,
              middle); // sort first section
    QuickSort(shorttype, theCollectionVector, middle + 1,
              bottom); // sort second section
  }
  return;
}

// Function to determine the partitions
// partitions the array and returns the middle subscript
int KM3SD::partition_CathodId(std::vector<KM3Hit *> *theCollectionVector,
                                int top, int bottom) {
  int x = (*theCollectionVector)[top]->GetCathodId();
  int i = top - 1;
  int j = bottom + 1;
  KM3Hit *temp;
  do {
    do {
      j--;
    } while (x < (*theCollectionVector)[j]->GetCathodId());

    do {
      i++;
    } while (x > (*theCollectionVector)[i]->GetCathodId());

    if (i < j) {
      temp = (*theCollectionVector)[i];
      (*theCollectionVector)[i] = (*theCollectionVector)[j];
      (*theCollectionVector)[j] = temp;
    }
  } while (i < j);
  return j; // returns middle subscript
}

// Function to determine the partitions
// partitions the array and returns the middle subscript
int KM3SD::partition_Time(std::vector<KM3Hit *> *theCollectionVector,
                            int top, int bottom) {
  double x = (*theCollectionVector)[top]->GetTime();
  int i = top - 1;
  int j = bottom + 1;
  KM3Hit *temp;
  do {
    do {
      j--;
    } while (x < (*theCollectionVector)[j]->GetTime());

    do {
      i++;
    } while (x > (*theCollectionVector)[i]->GetTime());

    if (i < j) {
      temp = (*theCollectionVector)[i];
      (*theCollectionVector)[i] = (*theCollectionVector)[j];
      (*theCollectionVector)[j] = temp;
    }
  } while (i < j);
  return j; // returns middle subscript
}

// the angular acceptance is according to the MultiPMT OM (WPD Document January
// 2011)
// doing linear interpolation
// if shapespherical==true then it is from parametrization and take into account
// only the <<experimental>> angular acceptance
// if shapespherical==false then it is not from param and take also the
// simulated angular acceptance of the cathod shape
bool KM3SD::AcceptAngle(double cosangle, double CathodRadius,
                          double CathodHeight, bool shapespherical) {
  static G4MaterialPropertyVector *Ang_Acc = NULL;
  static double MinCos_Acc = -1.0;
  static double MaxCos_Acc = 0.25;
  if (Ang_Acc == NULL) {
    G4Material *aMaterial = G4Material::GetMaterial("Cathod");
    Ang_Acc = aMaterial->GetMaterialPropertiesTable()->GetProperty(
        "ANGULAR_ACCEPTANCE");
    MinCos_Acc = Ang_Acc->GetMinLowEdgeEnergy();
    MaxCos_Acc = Ang_Acc->GetMaxLowEdgeEnergy();
  }

  if (cosangle > MaxCos_Acc)
    return false;
  if (cosangle < MinCos_Acc)
    return true;
  double AngAcc = Ang_Acc->Value(cosangle);

  // the following is the simulated angular acceptance (from the shape of the
  // photocathod)
  // A cylinder of radius R and height d has angular acceptance vs costh as
  // a(x)=abs(x)+(2*d/(pi*R))*sqrt(1-x*x)
  double AngularAccSim;

  if (!shapespherical) {
    AngularAccSim = fabs(cosangle) +
                    (2.0 * CathodHeight / (pi * CathodRadius)) *
                        sqrt(1 - cosangle * cosangle);
    AngAcc /= AngularAccSim;
  }

  if (G4UniformRand() <= AngAcc)
    return true;
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......