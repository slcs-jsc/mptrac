# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# List of executables...
EXC = atm2grid atm_conv atm_dist atm_init atm_select atm_split atm_stat cape day2doy doy2day jsec2time met_conv met_lapse met_map met_prof met_sample met_spec met_subgrid met_zm time2jsec sedi tnat trac tropo tropo_clim tropo_sample tropo_zm wind

# List of tests...
TESTS = atm_test compress_test interoper_test met_test tools_test trac_test tropo_test

# Installation directory...
DESTDIR ?= ../bin

# Include directories...
INCDIR += -I ../libs/build/include

# Library directories...
LIBDIR += -L ../libs/build/lib -L ../lib/build/lib64

# Compiler (gcc, nvc, or clang)...
COMPILER ?= gcc

# Compile for MPI...
MPI ?= 0

# Compile for GPUs...
GPU ?= 0

# GPU memory pinning for GH200...
GPU_PIN ?= 0

# Static compilation...
STATIC ?= 1

# Optimization flags...
OPTI ?= -O3

# Optimization information...
INFO ?= 0

# Compile for profiling...
PROF ?= 0

# Compile for coverage testing...
COV ?= 0

# Enable NVTX markers...
NVTX ?= 0

# Enable cuRAND random number generators...
CURAND ?= 0

# Enable async IO for GPUs (experimental)...
ASYNCIO ?= 0

# Enable KPP chemistry module (experimental)...
KPP ?= 0

# Enable Thrust sorting code (experimental)...
THRUST ?= 0

# Enable zfp compression (experimental)...
ZFP ?= 0

# Enable zstd compression (experimental)...
ZSTD ?= 0

# Enable cmultiscale compression (experimental)...
CMS ?= 0

# Path to cmultiscale library (experimental)...
CMSDIR ?= /home/lars/wrk/cmultiscale/cmultiscale

# -----------------------------------------------------------------------------
# Set compiler flags...
# -----------------------------------------------------------------------------

# Get git version...
VERSION ?= $(shell git describe --abbrev=6 --dirty --always --tags)

# Set CC and CFLAGS...
CFLAGS += $(INCDIR) $(DEFINES) -DVERSION=\"$(VERSION)\" $(OPTI) -g
ifeq ($(COMPILER),gcc)
  CC=gcc
  CFLAGS += -DHAVE_INLINE -pedantic -Werror -Wall -W -Wmissing-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp
else ifeq ($(COMPILER),nvc)
  CC=nvc
  CFLAGS += -mp
else ifeq ($(COMPILER),clang)
  CC=clang
  CFLAGS += -DHAVE_INLINE -pedantic -Werror -Wall -W -Wmissing-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wnested-externs -Wno-long-long -Wmissing-declarations -Wredundant-decls -Winline -fno-common -fshort-enums -fopenmp
else
  $(error Compiler $(COMPILER) not supported)
endif

# Compile for GPUs...
ifeq ($(GPU),1)
  ifeq ($(COMPILER),gcc)
    CFLAGS += -fopenacc -foffload=-lm
  else ifeq ($(COMPILER),nvc)
    CFLAGS += -acc
    ifeq ($(GPU_PIN),1)
      CFLAGS += -gpu=pinned
    endif
  else ifeq ($(COMPILER),clang)
    $(error OpenACC is not supported by the clang compiler)
  endif
endif

# Set LDFLAGS...
LDFLAGS += $(LIBDIR) -lgsl -lgslcblas -lnetcdf -lm
ifeq ($(STATIC),1)
  LDFLAGS += -lhdf5_hl -lhdf5 -lcurl -lsz -laec -lz -ldl -lm -lzstd
endif
ifeq ($(GPU),1)
  LDFLAGS += -L $(CUDA_PATH)/lib64 -lcudart
endif

# Compile with MPI...
ifeq ($(MPI),1)
  CC = mpicc
  CFLAGS += -DMPI
endif

# Optimization information...
ifeq ($(INFO),1)
  ifeq ($(COMPILER),gcc)
    CFLAGS += -fopt-info
  else ifeq ($(COMPILER),nvc)
    CFLAGS += -Minfo=all
  else ifeq ($(COMPILER),clang)
    $(error INFO flag not supported by clang compiler)
  endif
endif

# Compile for profiling...
ifeq ($(PROF),1)
  CFLAGS += -pg
endif

# Compile for coverage testing...
ifeq ($(COV),1)
  CFLAGS += --coverage
endif

# NVTX...
ifeq ($(NVTX),1)
  CFLAGS  += -DNVTX
  LDFLAGS += -lnvToolsExt
endif

# Compile with cuRAND...
ifeq ($(CURAND),1)
  ifeq ($(GPU),1)
    CFLAGS += -DCURAND
    LDFLAGS += -lcurand
  else
    $(error cuRAND requires compilation for GPUs)
  endif
endif

# asyncIO...
ifeq ($(ASYNCIO),1)
  ifneq ($(GPU),1)
    $(error Async-I/O requires compilation for GPUs)
  endif
  CFLAGS += -DASYNCIO
endif

# KPP chemistry...
ifeq ($(KPP),1)
  CFLAGS += -DKPP
  LDFLAGS += -lkpp
endif

# Thrust... 
ifeq ($(THRUST),1)
  ifeq ($(GPU),1)
    CFLAGS += -DTHRUST -c++libs -acc=gpu -gpu=cc80 -fast
    LDFLAGS += -lthrustsort_gpu
  else
    CFLAGS += -DTHRUST 
    LDFLAGS += -lthrustsort_cpu -L$(CUDA_PATH)/lib64 -lcudart -lstdc++ -lgomp 
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

# cmultiscale compression...
ifeq ($(CMS),1)
  LIBDIR += -L $(CMSDIR)/builddir/
  LDFLAGS += -lcmultiscale /usr/lib/x86_64-linux-gnu/libcurl.so.4
  ifeq ($(STATIC),1)
    LDFLAGS += -lstdc++
  endif
  INCDIR += -I $(CMSDIR)/src/cms_module/
  CFLAGS += -DCMS
endif

# Static compilation...
ifeq ($(STATIC),1)
  ifeq ($(GPU),1)
    $(error Static compilation does not work for GPU)
  endif
  ifeq ($(THRUST),1)
    $(error Static compilation does not work for THRUST)
  endif
  ifeq ($(MPI),1)
    $(error Static compilation does not work for MPI)
  endif
  ifeq ($(COMPILER),clang)
    $(error Static compilation does not work for clang compiler)
  endif
  CFLAGS += -static
endif

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

.PHONY : all check clean coverage cppcheck dist doxygen indent install lib lizard mkdocs strip uninstall

all: $(EXC)
	rm -f *~

$(EXC): %: %.c mptrac.o
	$(CC) $(CFLAGS) -o $@ $< mptrac.o $(LDFLAGS)

mptrac.o: mptrac.c mptrac.h Makefile
	$(CC) $(CFLAGS) -c -o mptrac.o mptrac.c

check: $(TESTS)

$(TESTS) gpu_test kpp_test: all
	@(echo "\n===== Running \"$@\" ... =====") ; \
	  cd ../tests/$@ ; ./run.sh \
	  && (echo "\n===== Test \"$@\" passed! =====") \
	  || (echo "\n===== Test \"$@\" failed! =====" ; exit 1)

clean:
	rm -rf $(EXC) lib* *.o *.gcda *.gcno *.gcov coverage* *~

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

doxygen:
	cd ../docs && doxygen

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

install:
	mkdir -p $(DESTDIR) && cp $(EXC) $(DESTDIR)

lib: mptrac.c mptrac.h Makefile
ifeq ($(STATIC),1)
	$(CC) $(CFLAGS) -c -o mptrac.o mptrac.c
	ar rcs libmptrac.a mptrac.o
else
	$(CC) $(CFLAGS) -c -fPIC -o mptrac_pic.o mptrac.c
	$(CC) mptrac_pic.o -shared -o libmptrac.so
endif

lizard:
	lizard -s cyclomatic_complexity

mkdocs:
	cd ../docs && mkdocs build

strip:
	strip $(EXC)

uninstall:
	cd $(DESTDIR) && rm $(EXC)
