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
  AnaTreeWriter.h

lib_LTLIBRARIES = \
   libanatreewriter.la

libanatreewriter_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `fastjet-config --libs`

libanatreewriter_la_SOURCES = \
  AnaTreeWriter.cc
  
libanatreewriter_la_LIBADD = \
  -lcalo_io \
  -lcalotrigger_io \
  -lcentrality_io \
  -ljetbase \
  -lfun4all \
  -lglobalvertex_io \
  -lglobalvertex \
  -lffarawobjects \
  -ljetbackground_io \
  -ljetbackground \
  -lphhepmc_io \
  -lphool \
  -lSubsysReco


BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals 

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libanatreewriter.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
