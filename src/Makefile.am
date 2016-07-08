common_header = devavxprng.h w256.hpp
noinst_PROGRAMS = dSFMTAVX2dc dSFMTAVX2eq
dSFMTAVX2dc_SOURCES = dSFMTAVX2dc.cpp dSFMTAVX2search.hpp \
dSFMTAVX2dc.h AlgorithmDSFMTEquidistribution.hpp \
AlgorithmCalcFixPoint.hpp Annihilate.hpp $(common_header)
dSFMTAVX2eq_SOURCES = dSFMTAVX2eq.cpp
AM_CXXFLAGS = -O2 -Wall -Wextra -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS
AM_LDFLAGS = -lMTToolBox
