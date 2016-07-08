common_header = devavxprng.h w256.hpp w512.hpp
noinst_DATA = Makefile.private dSFMTAVX2dc_mpi.cpp dSFMTAVX512Fdc_mpi.cpp
noinst_PROGRAMS = dSFMTAVX2dc dSFMTAVX2eq dSFMTAVX512Fdc dSFMTAVX512Feq
dSFMTAVX2dc_SOURCES = dSFMTAVX2dc.cpp dSFMTAVX2search.hpp \
dSFMTAVX2dc.h AlgorithmDSFMTEquidistribution.hpp \
AlgorithmCalcFixPoint.hpp Annihilate.hpp $(common_header)
dSFMTAVX2eq_SOURCES = dSFMTAVX2eq.cpp
dSFMTAVX512Fdc_SOURCES = dSFMTAVX512Fdc.cpp dSFMTAVX512Fsearch.hpp \
dSFMTAVX512Fdc.h $(common_header)
dSFMTAVX512Feq_SOURCES = dSFMTAVX512Feq.cpp
AM_CXXFLAGS = -O2 -Wall -Wextra -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS
AM_LDFLAGS = -lMTToolBox
