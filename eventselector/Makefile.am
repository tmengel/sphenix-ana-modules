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
  CentCut.h \
  EventCut.h \
  EventCutReport.h \
  EventSelector.h \
  LeadJetCut.h \
  LeadJetHook.h \
  MinBiasCut.h \
  TowerChi2Cut.h \
  ZVertexCut.h

lib_LTLIBRARIES = \
  libeventselection.la

libeventselection_la_SOURCES = \
  CentCut.cc \
  EventCutReport.cc \
  EventSelector.cc \
  LeadJetCut.cc \
  LeadJetHook.cc \
  MinBiasCut.cc \
  TowerChi2Cut.cc \
  ZVertexCut.cc

libeventselection_la_LIBADD = \
  -lphool \
  -ljetbase \
  -lg4dst \
  -lglobalvertex_io \
  -lglobalvertex \
  -lcentrality_io \
  -lSubsysReco

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libeventselection.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)