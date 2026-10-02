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

# the version number, in a single place (also read by
# .github/scripts/version.sh)
VERSION := $(shell cat VERSION 2>/dev/null)
ifeq ($(VERSION),)
  $(error cannot read the version number from ./VERSION)
endif

COMMON=-g -fno-exceptions -pthread

# Warnings of every recipe. The extra ones are known to every
# supported compiler (GCC 4.8.5 and later, clang), and swipe builds
# without any of them: once a warning of the DEBUG set below is fixed
# everywhere, its option moves here, so that it cannot come back
COMPILEOPT=-Wall -Wextra -Wpedantic -Wcast-align -Wcast-qual -Wconversion \
	-Wdouble-promotion -Wfloat-equal -Wformat=2 -Wnon-virtual-dtor \
	-Woverloaded-virtual -Wmissing-declarations -Wredundant-decls -Wshadow \
	-Wsign-conversion -Wswitch-default -Wold-style-cast -Wuninitialized \
	-Wunused -Wunused-macros -Wvla -Wzero-as-null-pointer-constant

# "make WERROR=1": warnings are errors (used by the CI; off by default,
# so that a compiler with new warnings can still build swipe)
ifdef WERROR
  COMPILEOPT += -Werror
endif

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
DEBUG_WARNINGS_COMMON=-Wdate-time -Wextra-semi -Wimplicit-fallthrough \
	-Wnull-dereference
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
ifeq ($(or $(RELEASE),$(DEBUG),$(PROFILE),$(COVERAGE)),)
  RELEASE := 1
endif

ifdef RELEASE
  # "make" or "make RELEASE=1": distributed binaries, benchmarks.
  # Reproducible: the debugging information names the source files
  # relative to the build directory, so that two builds of the same
  # sources in different directories are byte-identical
  # (.github/scripts/reproducible.sh)
  OPTIMIZATION=-O3 -DNDEBUG -fdebug-prefix-map=$(CURDIR)=.
else ifdef DEBUG
  # "make DEBUG=1": sanitizers and extended warnings (current GCC or
  # clang)
  OPTIMIZATION=-O0 -ggdb3 -D_GLIBCXX_DEBUG $(DEBUG_SANITIZER) $(DEBUG_WARNINGS)
  LINKOPT=$(DEBUG_SANITIZER)
else ifdef PROFILE
  # "make PROFILE=1": gprof (number of calls vs time per call)
  OPTIMIZATION=-pg -O1
  LINKOPT=-pg
else ifdef COVERAGE
  # "make COVERAGE=1": line and branch coverage with gcov
  OPTIMIZATION=-fprofile-arcs -ftest-coverage -O0
  LINKOPT=--coverage
  LIBS+=-lgcov
endif

# User variables (CXXFLAGS, CPPFLAGS, LDFLAGS, and LINKFLAGS, kept
# for compatibility) are appended after the flags above, so they can
# override them (e.g. make CXXFLAGS=-O2).
SWIPE_CXXFLAGS=$(STD) $(COMPILEOPT) $(COMMON) $(OPTIMIZATION) \
	$(VERSION_DEFINE) $(CPPFLAGS) $(CXXFLAGS)
SWIPE_LDFLAGS=$(COMMON) $(LINKOPT) $(LDFLAGS) $(LINKFLAGS)

PROG=swipe

all : swipe

# Installation directories (GNU conventions): make install PREFIX=...
# DESTDIR is prepended for staged installs (packaging)
PREFIX ?= /usr/local
exec_prefix := $(PREFIX)
datarootdir := $(PREFIX)/share
bindir := $(exec_prefix)/bin
mandir := $(datarootdir)/man
man1dir := $(mandir)/man1
bashcompdir ?= $(datarootdir)/bash-completion/completions
zshcompdir ?= $(datarootdir)/zsh/site-functions

MAN := man/swipe.1
BASH_COMPLETION := completion/swipe.bash
ZSH_COMPLETION := completion/_swipe

INSTALL ?= install
INSTALL_PROGRAM ?= $(INSTALL) -m 0755
INSTALL_DATA ?= $(INSTALL) -m 0644
MKDIR_P ?= $(INSTALL) -d

.PHONY : all clean distclean install install-completion uninstall

install : swipe $(MAN) install-completion
	$(MKDIR_P) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) swipe $(DESTDIR)$(bindir)/swipe
	$(MKDIR_P) $(DESTDIR)$(man1dir)
	$(INSTALL_DATA) $(MAN) $(DESTDIR)$(man1dir)/swipe.1

install-completion : $(BASH_COMPLETION) $(ZSH_COMPLETION)
	$(MKDIR_P) $(DESTDIR)$(bashcompdir)
	$(INSTALL_DATA) $(BASH_COMPLETION) $(DESTDIR)$(bashcompdir)/swipe
	$(MKDIR_P) $(DESTDIR)$(zshcompdir)
	$(INSTALL_DATA) $(ZSH_COMPLETION) $(DESTDIR)$(zshcompdir)/_swipe

uninstall :
	rm -f $(DESTDIR)$(bindir)/swipe
	rm -f $(DESTDIR)$(man1dir)/swipe.1
	rm -f $(DESTDIR)$(bashcompdir)/swipe
	rm -f $(DESTDIR)$(zshcompdir)/_swipe

clean :
	rm -f *.o *.d *~ $(PROG) gmon.out *.gcno *.gcda *.gcov

distclean : clean
	rm -f compile_commands.json

OBJS = options.o search_threads.o align_threads.o database.o asnparse.o align.o matrices.o \
	stats.o blastkar_partial.o hits.o query.o \
	search63.o search16.o search16s.o search7.o search7_ssse3.o search7_avx2.o

# Header dependencies are generated alongside each object (*.d
# files), so that editing any header rebuilds the right objects.
DEPFLAGS = -MMD -MP
DEPFILES = swipe.d $(OBJS:.o=.d)
-include $(DEPFILES)

DEPS = Makefile

# the version number (file VERSION) is only compiled into swipe.o,
# which defines swipe_name_and_version for the other files: a new
# version rebuilds swipe.o only
swipe.o : VERSION_DEFINE = -DSWIPE_VERSION='"$(VERSION)"'
swipe.o : VERSION

swipe : swipe.o $(OBJS)
	$(CXX) $(SWIPE_LDFLAGS) -o $@ $^ $(LIBS)

%.o : %.cc $(DEPS)
	$(CXX) $(SWIPE_CXXFLAGS) $(DEPFLAGS) -c -o $@ $<

search7_ssse3.o : search7.cc $(DEPS)
	$(CXX) -mssse3 $(SWIPE_CXXFLAGS) $(DEPFLAGS) -DSWIPE_SSSE3 -c -o $@ search7.cc

search7_avx2.o : search7_avx2.cc $(DEPS)
	$(CXX) -mavx2 $(SWIPE_CXXFLAGS) $(DEPFLAGS) -c -o $@ search7_avx2.cc
