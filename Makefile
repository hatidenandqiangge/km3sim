# $Id: GNUmakefile,v 1.1 1999/01/07 16:05:40 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := KM3Sim
G4TARGET := $(name)
G4EXLIB := true
G4WORKDIR := ./
CLHEP_BASE_DIR := /sps/km3net/users/tsirigot/HOURS/geant4.9.6.p02-install/
CLHEP_INCLUDE_DIR := /sps/km3net/users/tsirigot/HOURS/geant4.9.6.p02-install/include
CLHEP_LIB_DIR := /sps/km3net/users/tsirigot/HOURS/geant4.9.6.p02-install/lib64
CLHEP_LIB := G4clhep
G4SYSTEM := Linux-g++
#CXXFLAGS := -g -O0
CPPVERBOSE := true
#G4DEBUG := true

G4LIB_USE_EXPAT := true
G4LIB_BUILD_GDML := true
XERCESCROOT := /afs/cern.ch/sw/lcg/external/XercesC/3.1.1p2/x86_64-slc6-gcc47-opt
#in gdmllibs also load evt libs
GDMLLIBS := -levtio -lxerces-c -lpthread
G4LIB_BUILD_ZLIB := true
#CPPFLAGS += --static

CPPFLAGS += -I/sps/km3net/users/tsirigot/HOURS/v1r12seawiet/inc
LDFLAGS += -L/sps/km3net/users/tsirigot/HOURS/v1r12seawiet/lib -L/afs/cern.ch/sw/lcg/external/XercesC/3.1.1p2/x86_64-slc6-gcc47-opt/lib

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

ifndef G4INSTALL
G4INSTALL = /sps/km3net/users/tsirigot/HOURS/geant4.9.6.p02-install/share/Geant4-9.6.2/geant4make
endif

.PHONY: all
all: lib bin

G4INCLUDE = $(G4INSTALL)/../../../include/Geant4
G4LIB = $(G4INSTALL)/../../../lib64/Geant4-9.6.2


include $(G4INSTALL)/config/binmake.gmk

