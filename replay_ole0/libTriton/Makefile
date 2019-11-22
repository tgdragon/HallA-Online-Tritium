#------------------------------------------------------------------------------
# Names of source files and target libraries
# You do want to modify this section

# List all your source files here. They will be put into a shared library
# that can be loaded from a script.
# List only the implementation files (*.cxx). For every implementation file
# there must be a corresponding header file (*.h).
SRC  =	TriBCM.cxx
SRC +=	VETROCtdcModule.cxx TstF1TDCModule.cxx TdcDataEvtHandler.cxx
SRC +=	TriFadcCherenkov.cxx
SRC +=	TriFadcAerogelC.cxx
SRC +=	TriFadcRasteredBeam.cxx TriFadcUnRasteredBeam.cxx TriFadcRaster.cxx TriFadcBPM.cxx
SRC +=	TriFadcScin.cxx
SRC +=	TriFadcScinF1.cxx
SRC +=	TriFadcShower.cxx
SRC +=	TriFadcXscin.cxx
SRC +=	TriHRS.cxx TriSpectrometer.cxx TriSpectrometerDetector.cxx TriTrackingDetector.cxx \
	TriVDC.cxx TriVDCUVPlane.cxx TriVDCPlane.cxx TriVDCCluster.cxx TriVDCHit.cxx \
	TriVDCTimeToDistConv.cxx TriVDCUVTrack.cxx TriVDCTrackPair.cxx TriVDCTrackID.cxx \
	TriVDCWire.cxx TriTrack.cxx TriTrackInfo.cxx TriTrackingModule.cxx \
	TriNonTrackingDetector.cxx TriScintillator.cxx TriTrackID.cxx TriVDCAnalyticTTDConv.cxx \
	TriTriggerTime.cxx TriTrackProj.cxx TriXscin.cxx
SRC +=	TriScalerEvtHandler.cxx
SRC +=	V1495Module.cxx ClockCountEvtHandler.cxx
SRC +=	TriVDCeff.cxx
SRC +=	Tri_Beam_Eloss.cxx
SRC +=	Tri_ElossCorrection.cxx
SRC +=	Tri_Track_Eloss.cxx
SRC +=	Tritium_Beam_Eloss.cxx
SRC +=	Tritium_THaScaler100EvtHandler.cxx
SRC +=	Tritium_HRS.cxx
SRC +=	Tritium_TSScaler.cxx
SRC +=	Tritium_Track_Eloss.cxx
SRC +=	Tritium_Xscin.cxx

# Name of your package.
# The shared library that will be built will get the name lib$(PACKAGE).so
PACKAGE = Triton

# Name of the LinkDef file
LINKDEF = $(PACKAGE)_LinkDef.h

#------------------------------------------------------------------------------
# This part defines overall options and directory locations.
# Change as necessary,

# Compile with diagnostic code
export WITH_DEBUG = 1

# Compile debug version
#export DEBUG = 1

#------------------------------------------------------------------------------
# Directory locations. All we need to know is INCDIRS.
# INCDIRS lists the location(s) of the C++ Analyzer header (.h) files

# The following should work with both local installations and the
# Hall A counting house installation. For local installations, verify
# the setting of ANALYZER, or specify INCDIRS explicitly.

ifndef ANALYZER
  $(error $$ANALYZER environment variable not defined)
endif

INCDIRS  = $(wildcard $(addprefix $(ANALYZER)/, include))

#------------------------------------------------------------------------------
# Do not change anything  below here unless you know what you are doing

ifeq ($(strip $(INCDIRS)),)
  $(error No Analyzer header files found. Check $$ANALYZER)
endif

MACHINE := $(shell uname -s)
#ARCH    := linuxclang
ARCH    := linux
SOSUF   := so
ifeq ($(MACHINE),Darwin)
  ARCH := macosx
  SOSUF := dylib
endif

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTBIN      := $(shell root-config --bindir)
ROOTINC      := -I$(shell root-config --incdir)
CXX          := $(shell root-config --cxx)
CC           := $(shell root-config --cc)

INCLUDES      = $(addprefix -I, $(INCDIRS) ) -I$(shell pwd)

USERLIB       = lib$(PACKAGE).$(SOSUF)
USERDICT      = $(PACKAGE)Dict

LIBS          =
GLIBS         =

ifeq ($(ARCH),linux)
# Linux with egcs (>= RedHat 5.2)
ifdef DEBUG
  CXXFLAGS    = -g -O0
  LDFLAGS     = -g -O0
else
  CXXFLAGS    = -g -O
  LDFLAGS     = -g -O
  DEFINES     = -DNDEBUG
endif
CXXFLAGS     += -Wall -fPIC
LD            = g++
SOFLAGS       = -shared
DEFINES      += $(shell getconf LFS_CFLAGS)
endif

ifeq ($(ARCH),macosx)
# Mac OS X with gcc >= 3.x or clang++ >= 5
ifdef DEBUG
  CXXFLAGS   := -g -O0
  LDFLAGS    := -g -O0
  DEFINES    :=
else
  CXXFLAGS   := -O
  LDFLAGS    := -O
  DEFINES    := -DNDEBUG
endif
DEFINES      += -DMACVERS
CXXFLAGS     += -Wall -fPIC
LD           := $(CXX)
LDCONFIG     :=
SOFLAGS      := -shared -Wl,-undefined,dynamic_lookup
SONAME       := -Wl,-install_name,
CXXFLAGS     += -Wextra -Wno-missing-field-initializers -Wno-unused-parameter
endif

ifeq ($(ARCH),linuxclang)
# Linux with clang
CXX           = clang++
CC            = clang
ifdef DEBUG
  CXXFLAGS    = -g -O0
  LDFLAGS     = -g -O0
else
  CXXFLAGS    = -O
  LDFLAGS     = -O
endif
CXXFLAGS     += -Wall -Woverloaded-virtual -fPIC
LD            = g++
SOFLAGS       = -shared
endif

CXXFLAGS     += $(ROOTCFLAGS) $(INCLUDES) $(DEFINES)
LIBS         += $(ROOTLIBS) $(SYSLIBS)
GLIBS        += $(ROOTGLIBS) $(SYSLIBS)

ifdef WITH_DEBUG
DEFINES     += -DWITH_DEBUG
endif

ifndef PKG
PKG           = lib$(PACKAGE)
LOGMSG        = "$(PKG) source files"
else
LOGMSG        = "$(PKG) Software Development Kit"
endif
DISTFILE      = $(PKG).tar.gz

#------------------------------------------------------------------------------
OBJ           = $(SRC:.cxx=.o)
HDR           = $(SRC:.cxx=.h)
DEP           = $(SRC:.cxx=.d)
OBJS          = $(OBJ) $(USERDICT).o

all:		$(USERLIB)

$(USERLIB):	$(HDR) $(OBJS)
		$(LD) $(LDFLAGS) $(SOFLAGS) -o $@ $(OBJS)
		@echo "$@ done"

$(USERDICT).cxx: $(HDR) $(LINKDEF)
	@echo "Generating dictionary $(USERDICT)..."
	$(ROOTBIN)/rootcint -f $@ -c $(ROOTINC) $(INCLUDES) $(DEFINES) $^

install:	all
		$(error Please define install yourself)
# for example:
#		cp $(USERLIB) $(LIBDIR)

clean:
		rm -f *.o *~ $(USERLIB) $(USERDICT).* $(USERDICT)_rdict.pcm

realclean:	clean
		rm -f *.d

srcdist:
		rm -f $(DISTFILE)
		rm -rf $(PKG)
		mkdir $(PKG)
		cp -p $(SRC) $(HDR) $(LINKDEF) db*.dat README Makefile $(PKG)
		gtar czvf $(DISTFILE) --ignore-failed-read \
		 -V $(LOGMSG)" `date -I`" $(PKG)
		rm -rf $(PKG)

.PHONY: all clean realclean srcdist

.SUFFIXES:
.SUFFIXES: .c .cc .cpp .cxx .o .d

%.o:	%.cxx
ifeq ($(strip $(MAKEDEPEND)),)
	$(CXX) $(CXXFLAGS) -MMD -o $@ -c $<
	@mv -f $*.d $*.d.tmp
else
	$(CXX) $(CXXFLAGS) -o $@ -c $<
	$(MAKEDEPEND) $(ROOTINC) $(INCLUDES) $(DEFINES) -c $< > $*.d.tmp
endif
	@sed -e 's|.*:|$*.o:|' < $*.d.tmp > $*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

###

-include $(DEP)
