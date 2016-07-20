EXTRA_DIST = include/*
ACLOCAL_AMFLAGS = -I m4
SUBDIRS =
if MPI
  SUBDIRS += mpisrc
  EXTRA_DIST += src/*.cpp src/Makefile*
else
  SUBDIRS += src
  EXTRA_DIST += mpisrc/*.cpp mpisrc/Makefile*
endif
