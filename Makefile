# Makefile of pinc_MPI_bg

FC = mpif90

COMPILER = $(shell echo `mpif90 --version` | sed 's/ .*//')

ifeq ($(COMPILER), ifort)
  FCFLAGS=-O3 -real-size 64 -traceback -unroll=4 -ip
  MODULEFLAG=-module $(BUILD)
else
  FCFLAGS=-O3 -fdefault-real-8 -fbacktrace -funroll-loops
  MODULEFLAG= -J$(BUILD)
endif

LIPNAG =
LIBHYPRE = /home/atmodynamics/voelker/hypre/hypre-2.11.2/src/lib

# define directories for sources and binaries (GSV 072018)
BIN = ./bin
BUILD = ./build
SOURCE = ./src

OFILES =	types.o \
	mpi.o \
	timeScheme.o \
	algebra.o \
	atmosphere.o \
	init.o \
	debug.o \
	muscl.o \
	wkb.o \
	xweno.o \
	fluxes.o \
	boundary.o \
	poisson.o \
	update.o \
	output.o \
	finish.o \
	pinc.o \
	ice.o

# add build directory as prefix to path of %.o files
OBJ=$(addprefix $(BUILD)/, $(OFILES))

# general rules
$(BUILD)/%.o: $(SOURCE)/%.f90
	$(FC) $(FCFLAGS) $(MODULEFLAG) -c $< -o $@

# the main target
pinc:$(OBJ)
	$(FC) $(FCFLAGS) -o $(BIN)/pinc $(OBJ) $(LIPNAG77) -L$(LIBHYPRE) -lHYPRE

#-------------------
#  dependencies
#-------------------

# fluxes.f90
$(BUILD)/fluxes.o: $(BUILD)/types.o
$(BUILD)/fluxes.o: $(BUILD)/xweno.o
$(BUILD)/fluxes.o: $(BUILD)/muscl.o
$(BUILD)/fluxes.o: $(BUILD)/atmosphere.o
$(BUILD)/fluxes.o: $(BUILD)/algebra.o
$(BUILD)/fluxes.o: $(BUILD)/ice.o

# xweno.f90
$(BUILD)/xweno.o: $(BUILD)/types.o
$(BUILD)/xweno.o: $(BUILD)/debug.o

# debug.f90
$(BUILD)/debug.o: $(BUILD)/types.o
$(BUILD)/debug.o: $(BUILD)/atmosphere.o

# poisson.f90
$(BUILD)/poisson.o: $(BUILD)/types.o
$(BUILD)/poisson.o: $(BUILD)/mpi.o

$(BUILD)/mpi.o: $(BUILD)/fluxes.o
$(BUILD)/mpi.o: $(BUILD)/types.o
$(BUILD)/timeScheme.o: $(BUILD)/types.o
$(BUILD)/atmosphere.o: $(BUILD)/types.o
$(BUILD)/init.o: $(BUILD)/types.o
$(BUILD)/muscl.o: $(BUILD)/types.o
$(BUILD)/wkb.o: $(BUILD)/types.o
$(BUILD)/boundary.o: $(BUILD)/types.o
$(BUILD)/update.o: $(BUILD)/types.o
$(BUILD)/output.o: $(BUILD)/types.o
$(BUILD)/finish.o: $(BUILD)/types.o
$(BUILD)/ice.o: $(BUILD)/types.o

# test xweno_module
XOBJ = 	$(BUILD)/types.o $(BUILD)/xweno.o $(BUILD)/testXWENO.o $(BUILD)/debug.o
xweno:	$(XOBJ)
	$(FC) $(FCFLAGS) $(MODULEFLAG) -o testXWENO $(XOBJ)

# cleaning
TEMP = $(BUILD)/*.o $(BUILD)/*.mod $(BIN)/pinc
clean:
	rm -f $(TEMP)

##
#  OLD STUFF; DOESNT COMPILE DUE TO CHANGES IN THE CODE #

# program for file difference
DIFFOFILES = types.o \
	timeScheme.o \
	algebra.o \
	atmosphere.o \
	init.o \
	debug.o \
	muscl.o \
	wkb.o \
	xweno.o \
	fluxes.o \
	boundary.o \
	update.o \
	output.o \
	finish.o \
	algebra.o \
	l2diff.o

DIFFOBJ=$(addprefix $(BUILD)/, $(DIFFOFILES))

l2diff:$(DIFFOBJ)
	$(FC) $(FCFLAGS) -o l2diff $(DIFFOBJ)

# programme for making tec360 layout files
LAYOBJ = $(BUILD)/makeLayout.o $(BUILD)/types.o
layout: $(LAYOBJ)
	$(FC) $(FCFLAGS) $(MODULEFLAG) -o makeLayout $(LAYOBJ)

# test algebra_module
algebra: $(BUILD)/algebra.o
	$(FC) $(FCFLAGS) $(MODULEFLAG) -o testAlgebra $(SOURCE)/testAlgebra.f90 $(SOURCE)/algebra.f90


