INFO ?= 0
GPU ?= 0

# Set here the desired C compiler and its optimization options
ifeq ($(GPU),1)
  CC = nvc
	COPT = -O -mp -mcmodel=medium -acc -I ../../../libs/build/include -I ../../../src 
else
  CC   = gcc
	COPT = -O -Wall -fopenmp -I ../../../libs/build/include -I ../../../src
endif


# To create Matlab gateway routines
# Note: use $(CC) as the mex C compiler
MEX  = mex

HEADERS = chem_Global.h  chem_Parameters.h  chem_Sparse.h

SPSRC = chem_JacobianSP.c 

SPOBJ = chem_JacobianSP.o 

SRC =   chem_Main.c      chem_Integrator.c   \
	chem_Function.c  chem_Initialize.c   \
	chem_Jacobian.c  chem_LinearAlgebra.c\
	chem_Rates.c     chem_Hessian.c      \
	chem_Stoichiom.c chem_Util.c         \
	chem_Monitor.c

OBJ =  	chem_Main.o chem_Integrator.o chem_Function.o \
	chem_Jacobian.o chem_LinearAlgebra.o chem_Rates.o  \
	chem_Util.o  chem_JacobianSP.o \
	chem_Monitor.o

STOCHSRC = chem_Stochastic.c 
STOCHOBJ = chem_Stochastic.o 

all:    exe

exe:	$(HEADERS) $(SPOBJ) $(OBJ)
	$(CC) $(COPT) $(SPOBJ) $(OBJ) -lm -o chem.exe	

lib:  $(HEADERS) $(OBJ)

	ar rvs libkpp.a $(OBJ)

stochastic:$(HEADERS) $(SPOBJ) $(OBJ) $(STOCHOBJ)
	$(CC) $(COPT) $(SPOBJ) $(OBJ) $(STOCHOBJ) -lm \
	-o chem_stochastic.exe	

mex:    $(HEADERS) $(SPOBJ) $(OBJ)
	$(MEX) CC#$(CC) -O chem_mex_Fun.c     -lm $(SPOBJ) $(OBJ)
	$(MEX) CC#$(CC) -O chem_mex_Jac_SP.c  -lm $(SPOBJ) $(OBJ)
	$(MEX) CC#$(CC) -O chem_mex_Hessian.c -lm $(SPOBJ) $(OBJ)


clean:
	rm -f $(SPOBJ) $(OBJ) chem.exe chem_*.mexglx chem*.dat *.a

distclean:
	rm -f $(SPOBJ) $(OBJ) chem.exe chem*.dat \
	chem_*.c chem_*.h chem*.log chem_*.mexglx

chem_Monitor.o: chem_Monitor.c $(HEADERS)
	$(CC) $(COPT) -c $<

chem_JacobianSP.o: chem_JacobianSP.c $(HEADERS)
	$(CC) $(COPT) -c $<

chem_StoichiomSP.o: chem_StoichiomSP.c $(HEADERS)
	$(CC) $(COPT) -c $<

chem_Main.o: chem_Main.c chem_Initialize.o $(HEADERS) $(SPOBJ)
	$(CC) $(COPT) -c $<

chem_Integrator.o: chem_Integrator.c  $(HEADERS) $(SPOBJ)
	$(CC) $(COPT) -c $<

# chem_Initialize.o: chem_Initialize.c  $(HEADERS) $(SPOBJ)
# 	$(CC) $(COPT) -c $<

chem_Function.o: chem_Function.c  $(HEADERS) $(SPOBJ)
	$(CC) $(COPT) -c $<

chem_Stochastic.o: chem_Stochastic.c  $(HEADERS) $(SPOBJ)
	$(CC) $(COPT) -c $<

chem_Jacobian.o: chem_Jacobian.c  $(HEADERS) $(SPOBJ)
	$(CC) $(COPT) -c $<

chem_LinearAlgebra.o: chem_LinearAlgebra.c  $(HEADERS) $(SPOBJ)
	$(CC) $(COPT) -c $<

chem_Rates.o: chem_Rates.c  $(HEADERS) $(SPOBJ)
	$(CC) $(COPT) -c $<

chem_Hessian.o:  chem_Hessian.c $(HEADERS) $(SPOBJ)
	$(CC) $(COPT) -c $<

chem_Stoichiom.o: chem_Stoichiom.c  $(HEADERS) $(SPOBJ)
	$(CC) $(COPT) -c $<

chem_Util.o: chem_Util.c  $(HEADERS) $(SPOBJ)
	$(CC) $(COPT) -c $<

