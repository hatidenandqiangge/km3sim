# $Id: GNUmakefile,v 1.1 1999/01/07 16:05:40 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

CXX := mpic++

name := KM3Sim
G4TARGET := $(name)
G4EXLIB := true
G4WORKDIR := WORK
CLHEP_BASE_DIR := $(G4INSTALL)/
CLHEP_INCLUDE_DIR := $(G4INSTALL)/include
CLHEP_LIB_DIR := $(G4INSTALL)/lib64
CLHEP_LIB := G4clhep
G4SYSTEM := Linux-g++
#CXXFLAGS := -g -O0
CPPVERBOSE := true
#G4DEBUG := true

G4LIB_USE_EXPAT := true
G4LIB_BUILD_GDML := true
##XERCESCROOT := /afs/cern.ch/sw/lcg/external/XercesC/3.1.1p2/x86_64-slc6-gcc47-opt
#in gdmllibs also load evt libs
##GDMLLIBS := -levtio -lxerces-c -lpthread
GDMLLIBS := -lxerces-c -lpthread -lmpich -lmpichcxx

G4LIB_BUILD_ZLIB := true
#CPPFLAGS += --static

CPPFLAGS += -DG4MPI_SUPPORT
CPPFLAGS += -I/state/partition1/geant4.9.6.p02/source/examples/extended/parallel/MPI/source/include
CPPFLAGS += -I/usr/include/mpich-x86_64
LDFLAGS += -L/usr/lib64/mpich/lib

##LDFLAGS += -L/sps/km3net/users/tsirigot/HOURS/v1r12seawiet/lib -L/afs/cern.ch/sw/lcg/external/XercesC/3.1.1p2/x86_64-slc6-gcc47-opt/lib


CPPFLAGS += -DG4JUST_COUNT_PHOTONS
# Uncomment the following if you want to compile hadronic physics (includes muon photonuclear interaction).
# When hadronic physics is not included only basic particles are constructed
CPPFLAGS += -DG4HADRONIC_COMPILE

# Uncommend the following if you want to provide generation information for the optical photons
CPPFLAGS += -DG4TRACK_INFORMATION

# Uncomment the following if you want to have boundary proccess for optical photons
#CPPFLAGS += -DG4BOUNDARY_COMPILE

# comment out the following if you declare benthos and glass in geometrical construction of the OM
# else the absorption (only is considered in KM3SD)
# do not have this condition and glass or gell in geometrical construction
#CPPFLAGS += -DG4MY_TRANSPARENCIES

# uncomment the following if you are going to do energy, distance, photons parametrization for fit purposes
#CPPFLAGS += -DG4MYFIT_PARAMETERIZATION

# uncomment the following if you do not want to run with electromagnetic shower parametrization on
CPPFLAGS += -DG4DISABLE_PARAMETRIZATION

# uncomment the following if you want to activate Mie scattering proccess for optical photons
CPPFLAGS += -DG4ENABLE_MIE

# uncomment the following if you are going to do energy, distance, photons parametrization for em fast param purposes
#CPPFLAGS += -DG4MYEM_PARAMETERIZATION

# uncomment the following if you want Muon-Distance-Energy Parametrization
#CPPFLAGS += -DG4MYMUON_PARAMETERIZATION

# flag for creating Hadronic parametrization (photon tables)
#CPPFLAGS += -DG4MYHA_PARAMETERIZATION

# flag for creating Hadronic parametrization (muon tables)
# CPPFLAGS += -DG4MYHAMUONS_PARAMETERIZATION

# flag for printing header with cathod information
# CPPFLAGS += -DG4PRINT_HEADER


G4INSTALL := $(shell ls -d $(G4INSTALL)/share/G*)/geant4make

#ifndef G4INSTALL
#G4INSTALL = /sps/km3net/users/tsirigot/HOURS/geant4.9.6.p02-install/share/#Geant4-9.6.2/geant4make
#endif

.PHONY: all
all: lib bin

G4INCLUDE = $(G4INSTALL)/../../../include/Geant4
PREG4LIB = $(G4INSTALL)/../../../lib64/
G4LIB := $(shell ls -d $(PREG4LIB)/G*)
#G4LIB = $(G4INSTALL)/../../../lib64/Geant4-9.6.2

include $(G4INSTALL)/config/binmake.gmk

CXX :=/usr/lib64/mpich/bin/mpic++


