# Makefile of PincFloitMSGWAM
#
FC = mpif90
#FC = mpiifort   # recommended for Intel 2018 and onward

COMPILER = $(shell echo `$(FC) --version` | sed 's/ .*//')

ifeq ($(COMPILER), ifort)
 FCFLAGS=-O0 -g -check all -warn all -warn nounused -fpe0 -real-size 64 -traceback -unroll=4 -ip
#  FCFLAGS=-O3 -real-size 64 -traceback -unroll=4 -ip
 MODULEFLAG=-module $(BUILD)
else   # gcc
  #FCFLAGS=-O0 -g -fcheck=all -Wall -Wno-unused-variable -fdefault-real-8 -fbacktrace -funroll-loops -Wno-unused-dummy-argument -Wno-conversion-extra
 FCFLAGS=-O3 -fdefault-real-8 -fbacktrace -funroll-loops
 MODULEFLAG= -J$(BUILD)
endif

LIPNAG =
#LIBHYPRE ?= /home/atmodynamics/boeloeni/Hypre/hypre-2.11.2/src/lib
LIBHYPRE ?= /pf/b/b380792/Hypre/hypre-2.11.2/src/lib

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
	hypre_tools.o \
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

# pinc.f90
$(BUILD)/pinc.o: $(BUILD)/types.o
$(BUILD)/pinc.o: $(BUILD)/mpi.o
$(BUILD)/pinc.o: $(BUILD)/timeScheme.o
$(BUILD)/pinc.o: $(BUILD)/init.o
$(BUILD)/pinc.o: $(BUILD)/debug.o
$(BUILD)/pinc.o: $(BUILD)/wkb.o
$(BUILD)/pinc.o: $(BUILD)/output.o
$(BUILD)/pinc.o: $(BUILD)/xweno.o
$(BUILD)/pinc.o: $(BUILD)/atmosphere.o
$(BUILD)/pinc.o: $(BUILD)/boundary.o
$(BUILD)/pinc.o: $(BUILD)/fluxes.o
$(BUILD)/pinc.o: $(BUILD)/update.o
$(BUILD)/pinc.o: $(BUILD)/poisson.o
$(BUILD)/pinc.o: $(BUILD)/finish.o

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
$(BUILD)/poisson.o: $(BUILD)/hypre_tools.o
$(BUILD)/poisson.o: $(BUILD)/output.o

# wkb.f90
$(BUILD)/wkb.o: $(BUILD)/types.o
$(BUILD)/wkb.o: $(BUILD)/timeScheme.o
$(BUILD)/wkb.o: $(BUILD)/atmosphere.o
$(BUILD)/wkb.o: $(BUILD)/muscl.o

# mpi.f90
$(BUILD)/mpi.o: $(BUILD)/fluxes.o
$(BUILD)/mpi.o: $(BUILD)/types.o

# timeScheme.f90
$(BUILD)/timeScheme.o: $(BUILD)/types.o

# atmosphere.f90
$(BUILD)/atmosphere.o: $(BUILD)/types.o

# init.f90
$(BUILD)/init.o: $(BUILD)/types.o
$(BUILD)/init.o: $(BUILD)/ice.o
$(BUILD)/init.o: $(BUILD)/atmosphere.o
$(BUILD)/init.o: $(BUILD)/mpi.o
$(BUILD)/init.o: $(BUILD)/boundary.o

# muscl.f90
$(BUILD)/muscl.o: $(BUILD)/types.o

# boundary.f90
$(BUILD)/boundary.o: $(BUILD)/types.o

# update.f90
$(BUILD)/update.o: $(BUILD)/types.o
$(BUILD)/update.o: $(BUILD)/poisson.o
$(BUILD)/update.o: $(BUILD)/boundary.o

# output.f90
$(BUILD)/output.o: $(BUILD)/types.o

# finish.f90
$(BUILD)/finish.o: $(BUILD)/types.o

# ice.f90
$(BUILD)/ice.o: $(BUILD)/types.o
$(BUILD)/ice.o: $(BUILD)/atmosphere.o


# cleaning
TEMP = $(BUILD)/*.o $(BUILD)/*.mod $(BIN)/pinc
clean:
	rm -f $(TEMP)

##
#  OLD STUFF; DOESNT COMPILE DUE TO CHANGES IN THE CODE #

#  All obsolete below

## program for file difference
#DIFFOFILES = types.o \
#	timeScheme.o \
#	algebra.o \
#	atmosphere.o \
#	init.o \
#	debug.o \
#	muscl.o \
#	wkb.o \
#	xweno.o \
#	fluxes.o \
#	boundary.o \
#	update.o \
#	output.o \
#	finish.o \
#	algebra.o \
#	l2diff.o
#
#DIFFOBJ=$(addprefix $(BUILD)/, $(DIFFOFILES))
#
#l2diff:$(DIFFOBJ)
#	$(FC) $(FCFLAGS) -o l2diff $(DIFFOBJ)
#
## programme for making tec360 layout files
#LAYOBJ = $(BUILD)/makeLayout.o $(BUILD)/types.o
#layout: $(LAYOBJ)
#	$(FC) $(FCFLAGS) $(MODULEFLAG) -o makeLayout $(LAYOBJ)
#
## test algebra_module
#algebra: $(BUILD)/algebra.o
#	$(FC) $(FCFLAGS) $(MODULEFLAG) -o testAlgebra $(SOURCE)/testAlgebra.f90 $(SOURCE)/algebra.f90
#
## test xweno_module
#XOBJ = 	types.o xweno.o testXWENO.o debug.o
#xweno:	$(XOBJ)
#	$(FC) $(FCFLAGS) -o testXWENO $(XOBJ)
#
