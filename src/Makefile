# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Executables...
EXC = center cluster dist extract init jsec2time match met_map met_prof met_sample met_zm smago split time2jsec trac wind

# Library directories...
LIBDIR = -L ../lib/build/lib

# Include directories...
INCDIR = -I ../lib/build/include

# MPI...
#MPI = 1

# -----------------------------------------------------------------------------
# Set flags for GNU compiler...
# -----------------------------------------------------------------------------

# Compiler...
ifdef MPI
  CC = mpicc
else
  CC = gcc
endif

# CFLAGS...
CFLAGS = $(INCDIR) -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED -pedantic -Werror -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp -O3 -g

# LDFLAGS...
LDFLAGS = $(LIBDIR) -lgsl -lgslcblas -lnetcdf -lm 

# MPI...
ifdef MPI
  CFLAGS += -DMPI
else
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

bak:
	mkdir -p ../bak && zip ../bak/mptrac_`date +"%y%m%d%H%M"`.zip Doxyfile Makefile *.c *.h

clean:
	rm -f $(EXC) *.o *~

doc:
	mkdir -p ../doc && doxygen && cd ../doc/latex && make && mv refman.pdf ..

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h
