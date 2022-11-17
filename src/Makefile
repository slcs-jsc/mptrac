# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# List of executables...
EXC = atm_conv atm_dist atm_init atm_select atm_split atm_stat day2doy doy2day jsec2time met_conv met_lapse met_map met_prof met_sample met_spec met_subgrid met_zm time2jsec sedi tnat trac tropo tropo_sample tropo_zm wind

# List of tests...
TESTS = atm_test compress_test met_test tools_test trac_test tropo_test

# Installation directory...
DESTDIR ?= ../bin

# Include directories...
INCDIR += -I ../libs/build/include

# Library directories...
LIBDIR += -L ../libs/build/lib -L ../lib/build/lib64

# Compile for GPUs...
GPU ?= 0

# Compile for MPI...
MPI ?= 0

# Static compilation...
STATIC ?= 1

# ASYNCIO enables asunc IO for GPUs

ASYNCIO ?= 0

# Optimization flags...
OPT ?= -O3

# Optimization information...
INFO ?= 0

# Compile for profiling...
PROF ?= 0

# Compile for coverage report...
COV ?= 0

# Enable NVTX markers...
NVTX ?= 0

# Enable Thrust library...
THRUST ?= 0

# Enable zfp compression...
ZFP ?= 0

# Enable zstd compression...
ZSTD ?= 0

# Enable UVW cache...
UVW ?= 1

# -----------------------------------------------------------------------------
# Set compiler flags...
# -----------------------------------------------------------------------------

# Get git version...
VERSION ?= $(shell git describe --abbrev=4 --dirty --always --tags)

# Set CC and CFLAGS...
CFLAGS += $(INCDIR) $(DEFINES) -DVERSION=\"$(VERSION)\" $(OPT) -g
ifeq ($(GPU),1)
  CC = nvc
  CFLAGS += -mp -mcmodel=medium -acc
else ifeq ($(CC),nvc)
  CFLAGS += -mp -mcmodel=medium
else
  CC = gcc
  CFLAGS += -DHAVE_INLINE -pedantic -Werror -Wall -W -Wmissing-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp
endif

# Set LDFLAGS...
LDFLAGS += $(LIBDIR) -lgsl -lgslcblas -lnetcdf -lm
ifeq ($(STATIC),1)
  LDFLAGS += -lhdf5_hl -lhdf5 -lm -ldl -lz
endif

# Additional flags for GPUs...
ifeq ($(GPU),1)
  LDFLAGS += -L $(CUDA_PATH)/lib64 -lcudart -lcurand
endif

# Compile with MPI...
ifeq ($(MPI),1)
  CC = mpicc
  CFLAGS += -DMPI
endif

# Optimization information...
ifeq ($(INFO),1)
  ifeq ($(CC),nvc)
    CFLAGS += -Minfo=all
  else
    CFLAGS += -fopt-info
  endif
endif

# Compile for profiling...
ifeq ($(PROF),1)
  CFLAGS += -pg
endif

# Compile for coverage...
ifeq ($(COV),1)
  CFLAGS += -ftest-coverage -fprofile-arcs
endif

# NVTX...
ifeq ($(NVTX),1)
  CFLAGS  += -DNVTX
  LDFLAGS += -lnvToolsExt
endif

# Thrust...
ifeq ($(THRUST),1)
  ifeq ($(GPU),1)
    CFLAGS += -DTHRUST -c++libs
    LDFLAGS += -lthrustsort
  else
    $(error Thrust library only available for GPU)
  endif
endif

# zfp compression...
ifeq ($(ZFP),1)
  LDFLAGS += -lzfp
  CFLAGS += -DZFP
endif

# zstd compression...
ifeq ($(ZSTD),1)
  LDFLAGS += -lzstd
  CFLAGS += -DZSTD
endif

# UVW cache...
ifeq ($(UVW),1)
  CFLAGS += -DUVW
endif

# Static compilation...
ifeq ($(STATIC),1)
  ifeq ($(GPU),1)
    $(error Static compilation does not work for GPU)
  endif
  ifeq ($(MPI),1)
    $(error Static compilation does not work for MPI)
  endif
  CFLAGS += -static
endif

# asyncIO
ifeq ($(ASYNCIO),1)
  CFLAGS += -DASYNCIO
endif

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

.PHONY : all check clean coverage cppcheck dist doc indent \
           install lizard strip uninstall

all: $(EXC)
	rm -f *~

$(EXC): %: %.c libtrac.o
	$(CC) $(CFLAGS) -o $@ $< libtrac.o $(LDFLAGS)

check: $(TESTS)

$(TESTS): all
	cd ../tests/$@ ; ./run.sh \
	  && (echo "\nTest \"$@\" passed!") \
	  || (echo "\nTest \"$@\" failed!" ; exit 1)

libtrac.o: libtrac.c libtrac.h Makefile
	$(CC) $(CFLAGS) -c -o libtrac.o libtrac.c

clean:
	rm -rf $(EXC) *.o *.gcda *.gcno *.gcov coverage* *~

coverage:
	lcov --capture --directory . --output-file=coverage.info ; \
	genhtml coverage.info --legend --frames --output-directory=coverage
	gcov -abcfu *.c

cppcheck:
	cppcheck --enable=all ./

dist:
	cd ..  && tar cvfz mptrac_$(VERSION).tgz bin/* \
          docs/Doxyfile docs/refman.pdf example/* \
          libs/build.sh libs/*bz2 src/* tests/*

doc:
	cd ../docs && doxygen && cd latex && make && cp refman.pdf ..

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

install:
	mkdir -p $(DESTDIR) && cp $(EXC) $(DESTDIR)

lizard:
	lizard -s cyclomatic_complexity

strip:
	strip $(EXC)

uninstall:
	cd $(DESTDIR) && rm $(EXC)
