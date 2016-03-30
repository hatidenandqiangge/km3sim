#!/bin/bash
G4RADIOACTIVEDATA=/sps/km3net/users/tsirigot/HOURS/geant4.9.6.p02-install/share/Geant4-9.6.2/data/RadioactiveDecay3.6
G4LEVELGAMMADATA=/sps/km3net/users/tsirigot/HOURS/geant4.9.6.p02-install/share/Geant4-9.6.2/data/PhotonEvaporation2.3
G4LEDATA=/sps/km3net/users/tsirigot/HOURS/geant4.9.6.p02-install/share/Geant4-9.6.2/data/G4EMLOW6.32
G4NEUTRONHPDATA=/sps/km3net/users/tsirigot/HOURS/geant4.9.6.p02-install/share/Geant4-9.6.2/data/G4NDL4.2
G4SAIDXSDATA=/sps/km3net/users/tsirigot/HOURS/geant4.9.6.p02-install/share/Geant4-9.6.2/data/G4SAIDDATA1.1
G4REALSURFACEDATA=/sps/km3net/users/tsirigot/HOURS/geant4.9.6.p02-install/share/Geant4-9.6.2/data/RealSurface1.0
G4NEUTRONXSDATA=/sps/km3net/users/tsirigot/HOURS/geant4.9.6.p02-install/share/Geant4-9.6.2/data/G4NEUTRONXS1.2
G4PIIDATA=/sps/km3net/users/tsirigot/HOURS/geant4.9.6.p02-install/share/Geant4-9.6.2/data/G4PII1.3

export G4RADIOACTIVEDATA
export G4LEVELGAMMADATA
export G4LEDATA
export G4NEUTRONHPDATA
export G4SAIDXSDATA
export G4REALSURFACEDATA
export G4NEUTRONXSDATA
export G4PIIDATA


LD_LIBRARY_PATH=/sps/km3net/users/tsirigot/HOURS/geant4.9.6.p02-install/lib64:/afs/cern.ch/sw/lcg/external/XercesC/3.1.1p2/x86_64-slc6-gcc47-opt/lib:/sps/km3net/users/tsirigot/HOURS/v1r12seawiet/lib
export LD_LIBRARY_PATH

infile=$1
outfile=$2
seed=$3
./KM3Sim novrml $seed $outfile geo_orca_115strings_20mhorizontal_18OMs_6.0mvertical INPUTParametersRun_3500_scatter_WPD3.3_p0.0075_1.0 null null ANTARES_EVT_FORMAT $infile 
#> /dev/null
