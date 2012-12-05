# SWIPE
# 
# Smith-Waterman database searches with Inter-sequence Parallel Execution
# 
# Copyright (C) 2008-2012 Torbjørn Rognes, University of Oslo, 
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
# Contact: Torbjørn Rognes <torognes@ifi.uio.no>, 
# Department of Informatics, University of Oslo, 
# PO Box 1080 Blindern, NO-0316 Oslo, Norway

# Makefile for SWIPE

MPI_COMPILE=`mpicxx --showme:compile`
MPI_LINK=`mpicxx --showme:link`

COMMON=-g
#COMMON=-pg -g

LIBS=-lpthread
LINKFLAGS=$(COMMON)

# Intel options
CXX=icpc
CXXFLAGS=-Wall -Wno-missing-declarations -fast -xSSE2 $(COMMON)

# GNU options
#CXX=g++
#CXXFLAGS=-Wall -O3 -march=core2 $(COMMON)

PROG=swipe mpiswipe

OBJS = database.o asnparse.o align.o matrices.o \
	search7.o search16s.o search16.o search63.o \
	stats.o hits.o query.o

DEPS = swipe.h Makefile

.SUFFIXES:.o .cc

%.o : %.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

all : $(PROG)

mpiswipe.o : swipe.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -DMPISWIPE $(MPI_COMPILE) -c -o $@ swipe.cc

swipe : swipe.o $(OBJS)
	$(CXX) $(LINKFLAGS) -o $@ swipe.o $(OBJS) $(LIBS)

mpiswipe : mpiswipe.o $(OBJS)
	$(CXX) $(LINKFLAGS) -o $@ mpiswipe.o $(OBJS) $(LIBS) $(MPI_LINK)

clean :
	rm -f *.o *~ $(PROG) gmon.out
