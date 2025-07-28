AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -I`root-config --incdir`

lib_LTLIBRARIES = \
   libunderlyingevent_io.la \
   libunderlyingevent.la 

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

libunderlyingevent_io_la_LIBADD = \
  -lphool

libunderlyingevent_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `fastjet-config --libs`

libunderlyingevent_la_LIBADD = \
  libunderlyingevent_io.la \
  -ljetbase \
  -lfun4all \
  -lcalo_io \
  -lcalotrigger_io \
  -lcentrality_io \
  -lglobalvertex_io \
  -lphool \
  -lSubsysReco

pkginclude_HEADERS = \
  UEDefs.h \
  CaloWindowTowerReco.h \
  CaloWindowMap.h \
  CaloWindowMapv1.h \
  RandomConeTowerReco.h \
  RandomCone.h \
  RandomConev1.h 

ROOTDICTS = \
  CaloWindowMap_Dict.cc \
  CaloWindowMapv1_Dict.cc \
  RandomCone_Dict.cc \
  RandomConev1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  CaloWindowMap_Dict_rdict.pcm \
  CaloWindowMapv1_Dict_rdict.pcm \
  RandomCone_Dict_rdict.pcm \
  RandomConev1_Dict_rdict.pcm

libunderlyingevent_io_la_SOURCES = \
  $(ROOTDICTS) \
  CaloWindowMapv1.cc \
  RandomConev1.cc 

libunderlyingevent_la_SOURCES = \
  UEDefs.cc \
  CaloWindowTowerReco.cc \
  RandomConeTowerReco.cc 

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests
BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_underlyingevent_io \
  testexternals_underlyingevent


testexternals_underlyingevent_io_SOURCES = testexternals.cc
testexternals_underlyingevent_io_LDADD = libunderlyingevent_io.la

testexternals_underlyingevent_SOURCES = testexternals.cc
testexternals_underlyingevent_LDADD = libunderlyingevent.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
  
  