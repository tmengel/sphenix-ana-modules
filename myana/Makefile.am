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
  cppBash.h \
  AnaUtils.h \
  PlotUtils.h \
  sphenix_style.h

lib_LTLIBRARIES = \
   libmyana.la

libmyana_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib 

libmyana_la_SOURCES = \
  cppBash.cc \
  AnaUtils.cc \
  PlotUtils.cc \
  sphenix_style.cc
  
libmyana_la_LIBADD = \
  -lfun4all \
  -lSubsysReco 

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals 

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libmyana.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
