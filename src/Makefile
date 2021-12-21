# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Load these modules before compiling and running the code in Juelich:
#
# JURECA and JUWELS clusters:
# - enable MPI and disable STATIC flag
# - export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<...>/libs/build/lib
# ml purge
# ml GCC ParaStationMPI
#
# JURECA and JUWELS GPUs:
# - login to juwels-booster.fz-juelich.de
# - enable MPI and GPU flags and disable STATIC flag
# - export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<...>/libs/build/lib
# ml purge
# ml NVHPC ParaStationMPI

# List of executables...
EXC = atm_conv atm_dist atm_init atm_select atm_split atm_stat day2doy doy2day jsec2time lapse met_map met_prof met_sample met_spec met_zm time2jsec sedi subgrid tnat trac tropo tropo_sample wind

# Installation directory...
DESTDIR ?= ../bin

# Include directories...
INCDIR += -I ../libs/build/include

# Library directories...
LIBDIR += -L ../libs/build/lib -L ../lib/build/lib64

# Enable netCDF-4...
NC4 ?= 1

# Compile for GPUs...
GPU ?= 0

# Use managed memory for GPUs...
GPUMEM ?= 0

# Compile for MPI...
MPI ?= 0

# Static compilation...
STATIC ?= 1

# Optimization flags...
OPT ?= -O3

# Optimization information...
INFO ?= 0

# Compile for profiling...
PROF ?= 0

# Enable NVTX markers...
NVTX ?= 0

# -----------------------------------------------------------------------------
# Set compiler flags...
# -----------------------------------------------------------------------------

# Get git version...
VERSION ?= $(shell git describe --abbrev=4 --dirty --always --tags)

# Set CC and CFLAGS...
CFLAGS += $(INCDIR) $(DEFINES) -DVERSION=\"$(VERSION)\" $(OPT) -g

# Set CC and CFLAGS...
CFLAGS += $(INCDIR) $(DEFINES) $(OPT) -g
ifeq ($(GPU),1)
  CC = nvc
  CFLAGS += -mp -mcmodel=medium -acc
else ifeq ($(CC),nvc)
  CFLAGS += -mp -mcmodel=medium
else
  CC = gcc
  CFLAGS += -DHAVE_INLINE -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp
endif

# Set LDFLAGS...
LDFLAGS += $(LIBDIR) -lgsl -lgslcblas -lnetcdf -lm
ifeq ($(NC4),1)
  ifeq ($(STATIC),1)
    LDFLAGS += -lhdf5_hl -lhdf5 -lm -ldl -lz
  endif
endif

# Additional flags for GPUs...
ifeq ($(GPU),1)
  ifeq ($(GPUMEM),1)
    CFLAGS += -gpu=managed
  endif
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

# NVTX...
ifeq ($(NVTX),1)
  CFLAGS  += -DNVTX
  LDFLAGS += -lnvToolsExt
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

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

all: $(EXC)
	rm -f *~

$(EXC): %: %.c libtrac.o
	$(CC) $(CFLAGS) -o $@ $< libtrac.o $(LDFLAGS)

libtrac.o: libtrac.c libtrac.h Makefile
	$(CC) $(CFLAGS) -c -o libtrac.o libtrac.c

check: all
	cd ../example && ./run.sh

clean:
	rm -f $(EXC) *.o *~

dist:
	zip mptrac_`date +"%y%m%d%H%M"`.zip Doxyfile Makefile *.c *.h

doc:
	cd ../docs && doxygen && cd ../docs/latex && make

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

install: all
	mkdir -p $(DESTDIR) && cp $(EXC) $(DESTDIR)

strip:
	strip $(EXC)
