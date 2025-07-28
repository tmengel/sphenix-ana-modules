AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  CaloSpy.h \
  CaloTowerManip.h

lib_LTLIBRARIES = \
   libcalomanip.la

libcalomanip_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `fastjet-config --libs`

libcalomanip_la_SOURCES = \
  CaloSpy.cc \
  CaloTowerManip.cc 
  
libcalomanip_la_LIBADD = \
  -lcalo_io \
  -lcalotrigger_io \
  -lcentrality_io \
  -ljetbase \
  -ljetbackground_io \
  -ljetbackground \
  -lfun4all \
  -lglobalvertex_io \
  -lglobalvertex \
  -lphhepmc_io \
  -lphool \
  -lSubsysReco


BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals 

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libcalomanip.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
