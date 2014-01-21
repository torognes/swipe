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

COMMON=-g
#COMMON=-pg -g

COMPILEOPT=-Wall

LIBS=-lpthread

# Intel options
#CXX=icpc
#CXXFLAGS=$(COMPILEOPT) $(COMMON) -Wno-missing-declarations -fast
#LINKFLAGS=$(COMMON)

# GNU options
CXX=g++
CXXFLAGS=$(COMPILEOPT) $(COMMON) -O3
LINKFLAGS=$(COMMON)

PROG=swipe mpiswipe

all : $(PROG)

clean :
	rm -f *.o *~ $(PROG) gmon.out

OBJS = database.o asnparse.o align.o matrices.o \
	stats.o hits.o query.o \
	search63.o search16.o search16s.o search7.o search7_ssse3.o

DEPS = swipe.h Makefile

.SUFFIXES:.o .cc

swipe : swipe.o $(OBJS)
	$(CXX) $(LINKFLAGS) -o $@ $^ $(LIBS)

mpiswipe : mpiswipe.o $(OBJS)
	$(CXX) $(LINKFLAGS) -o $@ $^ $(LIBS) $(MPI_LINK)

%.o : %.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

mpiswipe.o : swipe.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -DMPISWIPE $(MPI_COMPILE) -c -o $@ swipe.cc

search7_ssse3.o : search7.cc $(DEPS)
	$(CXX) -mssse3 $(CXXFLAGS) -DSWIPE_SSSE3 -c -o $@ search7.cc
