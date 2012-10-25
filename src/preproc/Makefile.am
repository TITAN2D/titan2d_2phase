#
# automake template for preproc/ directory
#
# $Id: Makefile.am 12 2003-11-07 17:58:49Z kdalbey $
#

bin_PROGRAMS = hephaestus_materialnames hephaestus_preprocess

hephaestus_materialnames_SOURCES = materialnames.C

hephaestus_materialnames_CXXFLAGS = -I../gisapi/ -I../useful/

hephaestus_materialnames_LDADD = ../gisapi/libgisapi.a ../useful/libuseful.a 

hephaestus_preprocess_SOURCES = preprocess.C boundary.C createfunky.C element2.C \
	hilbert.C node.C \
	boundary.h element.h node.h \
	../header/FileFormat.h

hephaestus_preprocess_CXXFLAGS = -I../gisapi/ -I../useful/

hephaestus_preprocess_LDADD = ../gisapi/libgisapi.a ../useful/libuseful.a 


