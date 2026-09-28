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

COMMON=-g -pthread

# Warnings of every recipe. The extra ones are known to every
# supported compiler (GCC 4.8.5 and later, clang), and swipe builds
# without any of them: once a warning of the DEBUG set below is fixed
# everywhere, its option moves here, so that it cannot come back
COMPILEOPT=-Wall -Wextra -Wpedantic -Wcast-qual -Wdouble-promotion \
	-Wfloat-equal -Wformat=2 -Wnon-virtual-dtor -Woverloaded-virtual \
	-Wredundant-decls -Wshadow -Wswitch-default -Wuninitialized -Wunused \
	-Wunused-macros -Wvla -Wzero-as-null-pointer-constant

# language standard (swipe must build with GCC 4.8.5 and later)
STD=-std=c++11

LIBS=

# GNU options: g++, unless CXX is given (environment or command line)
ifeq ($(origin CXX),default)
  CXX=g++
endif

# Compiler identity: ask the preprocessor (works for g++, clang++,
# version-suffixed names, ccache wrappers and cross-compilers)
IS_CLANG := $(shell $(CXX) -x c++ -E -dM - < /dev/null 2>/dev/null | grep -c '__clang__')

# Extra warnings of the DEBUG recipe (current GCC and clang only: the
# oldest supported GCC rejects some of these options)
DEBUG_WARNINGS_COMMON=-Wcast-align -Wconversion -Wdate-time \
	-Wextra-semi -Wimplicit-fallthrough -Wnull-dereference \
	-Wold-style-cast -Wsign-conversion
DEBUG_WARNINGS_GCC=-Wduplicated-branches -Wduplicated-cond \
	-Wformat-overflow -Wlogical-op -Wuseless-cast
DEBUG_WARNINGS_CLANG=-Wcomma -Wassign-enum -Wover-aligned
ifneq ($(IS_CLANG),0)
  DEBUG_WARNINGS=$(DEBUG_WARNINGS_COMMON) $(DEBUG_WARNINGS_CLANG)
else
  DEBUG_WARNINGS=$(DEBUG_WARNINGS_COMMON) $(DEBUG_WARNINGS_GCC)
endif
DEBUG_SANITIZER=-fsanitize=undefined,address -fno-omit-frame-pointer

# Build recipes: exactly one is active, RELEASE by default. Objects
# are not tagged with the recipe: run "make clean" when switching.
ifeq ($(or $(RELEASE),$(DEBUG),$(PROFILE),$(COVERAGE),$(TRACE)),)
  RELEASE := 1
endif

ifdef RELEASE
  # "make" or "make RELEASE=1": distributed binaries, benchmarks
  OPTIMIZATION=-O3 -DNDEBUG
else ifdef DEBUG
  # "make DEBUG=1": sanitizers and extended warnings (current GCC or
  # clang). DEBUG is not defined: it enables the trace blocks (TRACE=1)
  OPTIMIZATION=-O0 -ggdb3 -D_GLIBCXX_DEBUG $(DEBUG_SANITIZER) $(DEBUG_WARNINGS)
  LINKOPT=$(DEBUG_SANITIZER)
else ifdef PROFILE
  # "make PROFILE=1": gprof (number of calls vs time per call)
  OPTIMIZATION=-pg -O1
  LINKOPT=-pg
else ifdef COVERAGE
  # "make COVERAGE=1": line and branch coverage with gcov
  OPTIMIZATION=-DCOVERAGE -fprofile-arcs -ftest-coverage -O0
  LINKOPT=--coverage
  LIBS+=-lgcov
else ifdef TRACE
  # "make TRACE=1": the #ifdef DEBUG trace blocks, printed to the output
  OPTIMIZATION=-O0 -ggdb3 -DDEBUG
endif

# User variables (CXXFLAGS, CPPFLAGS, LDFLAGS, and LINKFLAGS, kept
# for compatibility) are appended after the flags above, so they can
# override them (e.g. make CXXFLAGS=-O2).
SWIPE_CXXFLAGS=$(STD) $(COMPILEOPT) $(COMMON) $(OPTIMIZATION) $(CPPFLAGS) $(CXXFLAGS)
SWIPE_LDFLAGS=$(COMMON) $(LINKOPT) $(LDFLAGS) $(LINKFLAGS)

PROG=swipe

all : swipe

# Installation directories (GNU conventions): make install PREFIX=...
# DESTDIR is prepended for staged installs (packaging)
PREFIX ?= /usr/local
exec_prefix := $(PREFIX)
bindir := $(exec_prefix)/bin

INSTALL ?= install
INSTALL_PROGRAM ?= $(INSTALL) -m 0755
MKDIR_P ?= $(INSTALL) -d

.PHONY : all clean distclean install uninstall

install : swipe
	$(MKDIR_P) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) swipe $(DESTDIR)$(bindir)/swipe

uninstall :
	rm -f $(DESTDIR)$(bindir)/swipe

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
DEPFILES = swipe.d $(OBJS:.o=.d)
-include $(DEPFILES)

DEPS = Makefile

swipe : swipe.o $(OBJS)
	$(CXX) $(SWIPE_LDFLAGS) -o $@ $^ $(LIBS)

%.o : %.cc $(DEPS)
	$(CXX) $(SWIPE_CXXFLAGS) $(DEPFLAGS) -c -o $@ $<

search7_ssse3.o : search7.cc $(DEPS)
	$(CXX) -mssse3 $(SWIPE_CXXFLAGS) $(DEPFLAGS) -DSWIPE_SSSE3 -c -o $@ search7.cc
