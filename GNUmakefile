# $Id: GNUmakefile,v 1.2 2000-10-19 12:22:10 stanaka Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := firstTest
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

CPPFLAGS+= -g ## compile with debug information
########################### ROOT ################################
ifdef ROOTSYS
ifndef G4UI_USE_ROOT
CPPFLAGS += -DG4ANALYSIS_USE_ROOT $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS  = $(shell $(ROOTSYS)/bin/root-config --libs) 
ROOTLIBS += $(shell $(ROOTSYS)/bin/root-config --nonew --glibs) -lMinuit -lHtml
ROOTLIBS := $(filter-out -lNew,$(ROOTLIBS))
ROOTLIBS := $(filter-out -lpthread,$(ROOTLIBS))
EXTRALIBS += $(shell $(ROOTSYS)/bin/root-config --libs)
#LDLIBS += $(shell $(ROOTSYS)/bin/root-config --libs) 
LDLIBS += $(ROOTLIBS)
endif
endif
