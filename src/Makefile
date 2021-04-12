# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Load these modules before compiling and running the code:
#
# JURECA and JUWELS clusters:
#
# - enable MPI and disable STATIC flag
# - export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<...>/libs/build/lib
# ml purge
# ml GCC ParaStationMPI
#
# JURECA and JUWELS GPUs:
#
# - login to juwels-booster.fz-juelich.de
# - enable MPI and GPU flags and disable STATIC flag
# - export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<...>/libs/build/lib
# ml purge
# ml NVHPC ParaStationMPI

# List of executables...
EXC = atm_conv atm_dist atm_init atm_select atm_split atm_stat day2doy doy2day jsec2time lapse met_map met_prof met_sample met_spec met_zm time2jsec trac tropo tropo_sample

# Library directories...
LIBDIR += -L ../libs/build/lib -L ../lib/build/lib64

# Include directories...
INCDIR += -I ../libs/build/include

# Compile for GPUs...
GPU ?= 0

# Compile for MPI...
MPI ?= 0

# Compile for profiling...
PROF ?= 0

# Enable netCDF-4...
NC4 ?= 1

# Enable NVTX markers...
NVTX ?= 0

# Static compilation...
STATIC ?= 1

# -----------------------------------------------------------------------------
# Set compiler flags...
# -----------------------------------------------------------------------------

# Set CFLAGS...
CFLAGS += $(INCDIR) $(DEFINES) -g
ifeq ($(CC),pgcc)
  CFLAGS += -mp -mcmodel=medium -Minfo=all
else
  CFLAGS += -DHAVE_INLINE -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp
endif

# Set LDFLAGS...
LDFLAGS += $(LIBDIR) -lgsl -lgslcblas -lnetcdf -lm
ifeq ($(NC4),1)
  ifeq ($(STATIC),1)
    LDFLAGS += -lhdf5_hl -lhdf5 -lm -ldl -lz
  endif
endif

# Compile for GPUs...
ifeq ($(GPU),1)
  ifeq ($(CC),pgcc)
    CFLAGS += -acc -gpu=managed
  else
#    CFLAGS += -fopenacc -foffload=-lm
    CFLAGS += -fopenacc
  endif
  LDFLAGS += -L $(CUDA_PATH)/lib64 -lcudart -lcurand
endif

# Compile for MPI...
ifeq ($(MPI),1)
  CC = mpicc
  CFLAGS += -DMPI
endif

# Compile for profiling...
ifeq ($(PROF),1)
  CFLAGS += -O2 -pg
else
  CFLAGS += -O3
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

clean:
	rm -f $(EXC) *.o *~

dist:
	zip mptrac_`date +"%y%m%d%H%M"`.zip Doxyfile Makefile *.c *.h

doc:
	cd ../docs && doxygen && cd ../docs/latex && make

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

strip:
	strip $(EXC)
