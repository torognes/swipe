# SWIPE
# 
# Smith-Waterman database searches with Inter-sequence Parallel Execution
# 
# Copyright (C) 2008-2013 Torbjorn Rognes, University of Oslo, 
# Oslo University Hospital and Sencel Bioinformatics AS
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# Contact: Torbjorn Rognes <torognes@ifi.uio.no>, 
# Department of Informatics, University of Oslo, 
# PO Box 1080 Blindern, NO-0316 Oslo, Norway

# Makefile for SWIPE

MPI_COMPILE=`mpicxx --showme:compile`
MPI_LINK=`mpicxx --showme:link`

COMMON=-g -pthread
#COMMON=-pg -g

COMPILEOPT=-Wall -Wextra

# language standard (swipe must build with GCC 4.8.5 and later)
STD=-std=c++11

LIBS=

# GNU options: g++, unless CXX is given (environment or command line)
ifeq ($(origin CXX),default)
  CXX=g++
endif
OPTIMIZATION=-O3

# User variables (CXXFLAGS, CPPFLAGS, LDFLAGS, and LINKFLAGS, kept
# for compatibility) are appended after the flags above, so they can
# override them (e.g. make CXXFLAGS=-O2).
SWIPE_CXXFLAGS=$(STD) $(COMPILEOPT) $(COMMON) $(OPTIMIZATION) $(CPPFLAGS) $(CXXFLAGS)
SWIPE_LDFLAGS=$(COMMON) $(LDFLAGS) $(LINKFLAGS)

PROG=swipe mpiswipe

# mpiswipe (MPI version, needs mpicxx) is deprecated: it is no longer
# built by default, run "make mpiswipe" to build it
all : swipe

.PHONY : all clean distclean

clean :
	rm -f *.o *.d *~ $(PROG) gmon.out *.gcno *.gcda *.gcov

distclean : clean
	rm -f compile_commands.json

OBJS = database.o asnparse.o align.o matrices.o \
	stats.o hits.o query.o \
	search63.o search16.o search16s.o search7.o search7_ssse3.o

# Header dependencies are generated alongside each object (*.d
# files), so that editing any header, or blastkar_partial.c (included
# by stats.cc), rebuilds the right objects.
DEPFLAGS = -MMD -MP
DEPFILES = swipe.d mpiswipe.d $(OBJS:.o=.d)
-include $(DEPFILES)

DEPS = Makefile

swipe : swipe.o $(OBJS)
	$(CXX) $(SWIPE_LDFLAGS) -o $@ $^ $(LIBS)

mpiswipe : mpiswipe.o $(OBJS)
	$(CXX) $(SWIPE_LDFLAGS) -o $@ $^ $(LIBS) $(MPI_LINK)

%.o : %.cc $(DEPS)
	$(CXX) $(SWIPE_CXXFLAGS) $(DEPFLAGS) -c -o $@ $<

mpiswipe.o : swipe.cc $(DEPS)
	$(CXX) $(SWIPE_CXXFLAGS) $(DEPFLAGS) -DMPISWIPE $(MPI_COMPILE) -c -o $@ swipe.cc

search7_ssse3.o : search7.cc $(DEPS)
	$(CXX) -mssse3 $(SWIPE_CXXFLAGS) $(DEPFLAGS) -DSWIPE_SSSE3 -c -o $@ search7.cc
