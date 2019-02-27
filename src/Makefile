# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Load these modules before compiling and running the code:
#
# JURECA cluster:
# module load GCC/7.3.0 ParaStationMPI/5.2.1-1
#
# JURECA booster:
# module purge && module load Architecture/KNL GCC/8.2.0 ParaStationMPI/5.2.1-1

# List of executables...
EXC = atm_conv atm_dist atm_init atm_split atm_stat day2doy doy2day extract jsec2time met_map met_prof met_sample met_zm smago time2jsec trac wind

# Library directories...
LIBDIR = -L ../lib/build/lib -L ../lib/build/lib64

# Include directories...
INCDIR = -I ../lib/build/include

# Compile for KNL...
#KNL = 1

# Compile for MPI...
#MPI = 1

# Compile for profiling...
#PROF = 1

# -----------------------------------------------------------------------------
# Set flags for GNU compiler...
# -----------------------------------------------------------------------------

# Select compiler...
ifdef MPI
  CC = mpicc
else
  CC = gcc
endif

# Set CFLAGS...
CFLAGS = $(INCDIR) -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp -g

# Set LDFLAGS...
LDFLAGS = $(LIBDIR) -lgsl -lgslcblas -lnetcdf -lm

# Compile for KNL...
ifdef KNL
  CFLAGS += -march=knl -mtune=knl
endif

# Compile for MPI...
ifdef MPI
  CFLAGS += -DMPI
else
  CFLAGS += -static
endif

# Compile for profiling...
ifdef PROF
  CFLAGS += -O2 -pg
else
  CFLAGS += -O3
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
	mkdir -p ../doc && doxygen && cd ../doc/latex && make && mv refman.pdf ..

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

strip:
	strip $(EXC)
