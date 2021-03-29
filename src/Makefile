# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Load these modules before compiling and running the code:
#
# JUWELS cluster:
# - enable MPI and disable STATIC flag
# - export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<...>/libs/build/lib
# ml GCC ParaStationMPI
#
# JUWELS booster:
# - login to juwels-booster.fz-juelich.de
# - enable MPI and GPU flags and disable STATIC flag
# - export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<...>/libs/build/lib
# ml NVHPC ParaStationMPI
#
# JURECA-DC cluster:
# - see JUWELS cluster
#
# JURECA-DC GPUs:
# - see JUWELS booster
#
# JURECA-DC KNLs:
# ml purge
# ml Architecture/KNL GCC ParaStationMPI

# List of executables...
EXC = atm_conv atm_dist atm_init atm_select atm_split atm_stat day2doy doy2day jsec2time lapse met_map met_prof met_sample met_spec met_zm time2jsec trac tropo tropo_sample

# Library directories...
LIBDIR = -L ../libs/build/lib -L ../lib/build/lib64

# Include directories...
INCDIR = -I ../libs/build/include

# Compile for GPUs...
#GPU = 1

# Compile for KNL...
#KNL = 1

# Compile for MPI...
#MPI = 1

# Compile for profiling...
#PROF = 1

# Enable NVTX markers...
#USE_NVTX = 1

# Static compilation...
STATIC = 1

# -----------------------------------------------------------------------------
# Set flags for GNU compiler...
# -----------------------------------------------------------------------------

# Set compiler...
CC = gcc

# Set CFLAGS...
CFLAGS = $(INCDIR) -DHAVE_INLINE -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp -g

# Set LDFLAGS...
LDFLAGS = $(LIBDIR) -lgsl -lgslcblas -lnetcdf -lm

# Compile for GPUs...
ifdef GPU
  CC = pgcc
  CFLAGS = $(INCDIR) -mp -acc -gpu=managed -mcmodel=medium -Minfo=all -g
  LDFLAGS += -L $(CUDA_PATH)/lib64 -lcudart -lcurand
endif

# Compile for KNL...
ifdef KNL
  CFLAGS += -march=knl -mtune=knl
endif

# Compile for MPI...
ifdef MPI
  CC = mpicc
  CFLAGS += -DMPI
endif

# Compile for profiling...
ifdef PROF
  CFLAGS += -O2 -pg
else
  CFLAGS += -O3
endif

# NVTX...
ifdef USE_NVTX
  LDFLAGS += -lnvToolsExt
  CFLAGS  += -DUSE_NVTX
endif

# Static compilation...
ifdef STATIC
ifdef GPU
$(error Static compilation does not work for GPU)
endif
ifdef MPI
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
	mkdir -p ../docs && doxygen && cd ../docs/latex && make && mv refman.pdf ..

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

strip:
	strip $(EXC)
