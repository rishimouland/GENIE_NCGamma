#
# Makefile for GENIE Neutrino Generator
#
# Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL    = /bin/sh
NAME     = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

PACKAGE           = NuGamma
DICTIONARY        = _ROOT_DICT_NuGamma
LIBNAME           = libGNuGamma
EXTRA_EXT_LIBS    = cgaus.f cs.f90

all     : rootcint lib lib-link cs.o
install : install-inc install-lib

cs.o: cs.f90 diag.f90 cgaus.f
	$(FORT) -c cs.f90 cgaus.f -fPIC

# Include standard package makefile targets
#
include $(GENIE)/src/make/Make.std-package-targets

FORCE:

