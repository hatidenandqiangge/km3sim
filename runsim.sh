#!/bin/bash
G4RADIOACTIVEDATA=/afs/cern.ch/sw/lcg/external/geant4/9.5.p02/data/RadioactiveDecay3.4
G4LEVELGAMMADATA=/afs/cern.ch/sw/lcg/external/geant4/9.5.p02/data/PhotonEvaporation2.2
G4LEDATA=/afs/cern.ch/sw/lcg/external/geant4/9.5.p02/data/G4EMLOW6.23
G4NEUTRONHPDATA=/afs/cern.ch/sw/lcg/external/geant4/9.5.p02/data/G4NDL4.0
G4ABLADATA=/afs/cern.ch/sw/lcg/external/geant4/9.5.p02/data/G4ABLA3.0
G4REALSURFACEDATA=/afs/cern.ch/sw/lcg/external/geant4/9.5.p02/data/RealSurface1.0
G4NEUTRONXSDATA=/afs/cern.ch/sw/lcg/external/geant4/9.5.p02/data/G4NEUTRONXS1.1
G4PIIDATA=/afs/cern.ch/sw/lcg/external/geant4/9.5.p02/data/G4PII1.3

export G4RADIOACTIVEDATA
export G4LEVELGAMMADATA
export G4LEDATA
export G4NEUTRONHPDATA
export G4ABLADATA
export G4REALSURFACEDATA
export G4NEUTRONXSDATA
export G4PIIDATA

ntype=$1
runnumber=$2

LD_LIBRARY_PATH=/sps/km3net/users/tsirigot/HOURS/run1_ver4/WORK/tmp/Linux-g++/KM3Sim:/sps/km3net/users/tsirigot/HOURS/geant4.9.5.p02-install/lib64:/afs/cern.ch/sw/lcg/external/XercesC/3.1.1p2/x86_64-slc5-gcc43-opt/lib:/sps/km3net/users/tsirigot/HOURS/v1r12seawiet/lib
APOSTOLOS_PATH=/sps/km3net/users/tsirigot/ORCARunDense

echo ${APOSTOLOS_PATH}/bartol_genie_gene/run20/$ntype/bartol_NU_"$ntype"_Emin_1_Emax_30_seed_$(($runnumber))_run_20.evt.gz


if [ -f ${APOSTOLOS_PATH}/bartol_genie_gene/run20/$ntype/bartol_NU_"$ntype"_Emin_1_Emax_30_seed_$(($runnumber))_run_20.evt.gz ] ; then
    cp ${APOSTOLOS_PATH}/bartol_genie_gene/run20/$ntype/bartol_NU_"$ntype"_Emin_1_Emax_30_seed_$(($runnumber))_run_20.evt.gz ./
    gunzip bartol_NU_"$ntype"_Emin_1_Emax_30_seed_$(($runnumber))_run_20.evt.gz
  
    # Parsing is a wholly different story if pythia is involved, but I suppose
    # that I can just ignore that completely (for now)
    # 0      1      2          3                                                                    4                                                          5                                                     6    7    8                  9                                                       
    #        UI     seed       outfile                                                              geom                                                       param-file(physics)                                empara HAparam  antformat?      input-particles-nuflux                          
    ./KM3Sim novrml $runnumber simbartol_NU_"$ntype"_Emin_1_Emax_30_seed_$(($runnumber))_run_20.evt geo_orca_octagon_2181strings_3mhorizontal_51f_3.0mvertical INPUTParametersRun_Oscil_3500_scatter_WPD_p0.0075_1.0 null null ANTARES_EVT_FORMAT bartol_NU_"$ntype"_Emin_1_Emax_30_seed_$(($runnumber))_run_20.evt

#    gzip bartol_NU_"$ntype"_Emin_1_Emax_30_seed_$(($runnumber))_run_10.evt
else
    echo No such run for neutrino type $ntype 
fi

