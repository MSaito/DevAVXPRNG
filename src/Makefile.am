EXTRA_DIST = Makefile.private dSFMTAVX2dc_mpi.cpp dSFMTAVX512Fdc_mpi.cpp \
dSFMTAVX2eq_mpi.cpp dSFMTAVX512Feq_mpi.cpp \
SFMTAVX2dc_mpi.cpp SFMTAVX512Fdc_mpi.cpp \
SFMTAVX2eq_mpi.cpp SFMTAVX512Feq_mpi.cpp

mpi_exes = dSFMTAVX2dc_mpi dSFMTAVX512Fdc_mpi \
dSFMTAVX2eq_mpi dSFMTAVX512Feq_mpi \
SFMTAVX2dc_mpi SFMTAVX512Fdc_mpi \
SFMTAVX2eq_mpi SFMTAVX512Feq_mpi

common_header = devavxprng.h w256.hpp w512.hpp
noinst_PROGRAMS = dSFMTAVX2dc dSFMTAVX2eq dSFMTAVX512Fdc dSFMTAVX512Feq \
SFMTAVX2dc SFMTAVX2eq SFMTAVX512Fdc SFMTAVX512Feq
if MPI
  noinst_PROGRAMS += $(mpi_exes)
  dSFMTAVX2dc_mpi_SOURCES = dSFMTAVX2dc_mpi.cpp
  dSFMTAVX512Fdc_mpi_SOURCES = dSFMTAVX512Fdc_mpi.cpp
  dSFMTAVX2eq_mpi_SOURCES = dSFMTAVX2eq_mpi.cpp
  dSFMTAVX512Feq_mpi_SOURCES = dSFMTAVX512Feq_mpi.cpp
  SFMTAVX2dc_mpi_SOURCES = SFMTAVX2dc_mpi.cpp
  SFMTAVX512Fdc_mpi_SOURCES = SFMTAVX512Fdc_mpi.cpp
  SFMTAVX2eq_mpi_SOURCES = SFMTAVX2eq_mpi.cpp
  SFMTAVX512Feq_mpi_SOURCES = SFMTAVX512Feq_mpi.cpp
endif

dSFMTAVX2dc_SOURCES = dSFMTAVX2dc.cpp dSFMTAVX2search.hpp \
dSFMTAVX2dc.h AlgorithmDSFMTEquidistribution.hpp \
AlgorithmCalcFixPoint.hpp Annihilate.hpp $(common_header)
dSFMTAVX2eq_SOURCES = dSFMTAVX2eq.cpp
dSFMTAVX512Fdc_SOURCES = dSFMTAVX512Fdc.cpp dSFMTAVX512Fsearch.hpp \
dSFMTAVX512Fdc.h $(common_header)
dSFMTAVX512Feq_SOURCES = dSFMTAVX512Feq.cpp
SFMTAVX2dc_SOURCES = SFMTAVX2dc.cpp SFMTAVX2search.hpp \
SFMTAVX2dc.h AlgorithmSIMDEquidistribution.hpp $(common_header)
SFMTAVX2eq_SOURCES = SFMTAVX2eq.cpp
SFMTAVX512Fdc_SOURCES = SFMTAVX512Fdc.cpp SFMTAVX512Fsearch.hpp \
SFMTAVX512Fdc.h $(common_header)
SFMTAVX512Feq_SOURCES = SFMTAVX512Feq.cpp
AM_CXXFLAGS = -O2 -Wall -Wextra -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS
AM_LDFLAGS = -lMTToolBox
