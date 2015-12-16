# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = center dist extract init jsec2time met_map met_prof smago split time2jsec trac

# Library directories...
LIBDIR = -L../lib/build/lib -L../../project/lib/build/lib -L../../project/lib/build/lib64 -L/usr/local/gsl/1.14/lib -L/usr/local/netcdf/v4.3.0/lib

# Include directories...
INCDIR = -I../lib/build/include -I../../project/lib/build/include -I/usr/local/gsl/1.14/include -I/usr/local/netcdf/v4.3.0/include

# Linking...
STATIC = 1

# Profiling...
#PROF = 1

# MPI...
#MPI = 1

# -----------------------------------------------------------------------------
# Set flags...
# -----------------------------------------------------------------------------

# Compiler...
ifdef MPI
  CC = mpicc
else
  CC = gcc
endif

# CFLAGS...
CFLAGS = $(INCDIR) -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp

# LDFLAGS...
LDFLAGS = $(LIBDIR) -lgsl -lgslcblas -lnetcdf -lm 

# Profiling...
ifdef PROF
  CFLAGS += -O2 -g -pg 
else
  CFLAGS += -O3
endif

# Linking...
ifdef STATIC
  CFLAGS += -static
endif

# MPI...
ifdef MPI
  CFLAGS += -DMPI
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

bak:
	mkdir -p ../bak && zip ../bak/mptrac_`date +"%y%m%d%H%M"`.zip Doxyfile Makefile *.c *.h

clean:
	rm -f $(EXC) *.o *~

doc:
	mkdir -p ../doc && doxygen && cd ../doc/latex && make

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h
