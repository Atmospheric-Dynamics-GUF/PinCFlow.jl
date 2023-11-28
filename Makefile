# Makefile of PincFlow

# Set compiler.
FC = mpif90
# FC = mpiifort
COMPILER = $(shell echo `$(FC) --version` | sed 's/ .*//')

# Set flags.
ifeq ($(COMPILER), ifort)
  # FCFLAGS=-O0 -g -check all -real-size 64 -traceback -unroll=4 -ip
  FCFLAGS=-O3 -real-size 64 -traceback -unroll=4 -ip
  MODULEFLAG=-module $(BUILD)
else # GNU
  # FCFLAGS=-O0 -g -fcheck=all -Wall -Wno-unused-variable -fdefault-real-8 -fbacktrace -funroll-loops -Wno-unused-dummy-argument -Wno-conversion-extra
  FCFLAGS=-O3 -fdefault-real-8 -fbacktrace -funroll-loops -fallow-argument-mismatch -w
  MODULEFLAG= -J$(BUILD)
endif

# Define directories for sources and binaries.
BIN = $(shell mkdir -p ./bin) ./bin
BUILD = $(shell mkdir -p ./build) ./build
SOURCE = ./src

# Define *.o files.
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
	bicgstab_tools.o \
	poisson.o \
	update.o \
	output.o \
	finish.o \
	pinc.o \
	ice.o \
	sizeof.o \
	tracer.o \
	ice2.o \
	ice2_sub.o \
	optField.o

# Add build directory as prefix to path of *.o files.
OBJ=$(addprefix $(BUILD)/, $(OFILES))

# Set general rules.
$(BUILD)/%.o: $(SOURCE)/%.f90
	$(FC) $(FCFLAGS) $(MODULEFLAG) -c $< -o $@

# Set the main target.
pinc:$(OBJ)
	$(FC) $(FCFLAGS) -o $(BIN)/pinc $(OBJ)
	mkdir -p code
	cp -r $(SOURCE)/* ./code

# List dependencies.

# pinc.f90
$(BUILD)/pinc.o: $(BUILD)/types.o
$(BUILD)/pinc.o: $(BUILD)/sizeof.o
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
$(BUILD)/pinc.o: $(BUILD)/tracer.o
$(BUILD)/pinc.o: $(BUILD)/ice2.o

# fluxes.f90
$(BUILD)/fluxes.o: $(BUILD)/types.o
$(BUILD)/fluxes.o: $(BUILD)/xweno.o
$(BUILD)/fluxes.o: $(BUILD)/muscl.o
$(BUILD)/fluxes.o: $(BUILD)/atmosphere.o
$(BUILD)/fluxes.o: $(BUILD)/algebra.o
$(BUILD)/fluxes.o: $(BUILD)/ice.o
$(BUILD)/fluxes.o: $(BUILD)/sizeof.o

# xweno.f90
$(BUILD)/xweno.o: $(BUILD)/types.o
$(BUILD)/xweno.o: $(BUILD)/debug.o

# debug.f90
$(BUILD)/debug.o: $(BUILD)/types.o
$(BUILD)/debug.o: $(BUILD)/atmosphere.o

# poisson.f90
$(BUILD)/poisson.o: $(BUILD)/types.o
$(BUILD)/poisson.o: $(BUILD)/mpi.o
$(BUILD)/poisson.o: $(BUILD)/bicgstab_tools.o
$(BUILD)/poisson.o: $(BUILD)/output.o
$(BUILD)/poisson.o: $(BUILD)/sizeof.o

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
$(BUILD)/atmosphere.o: $(BUILD)/sizeof.o

# init.f90
$(BUILD)/init.o: $(BUILD)/types.o
$(BUILD)/init.o: $(BUILD)/ice.o
$(BUILD)/init.o: $(BUILD)/atmosphere.o
$(BUILD)/init.o: $(BUILD)/mpi.o
$(BUILD)/init.o: $(BUILD)/boundary.o
$(BUILD)/init.o: $(BUILD)/sizeof.o
$(BUILD)/init.o: $(BUILD)/tracer.o
$(BUILD)/init.o: $(BUILD)/ice2.o
$(BUILD)/init.o: $(BUILD)/ice2_sub.o
$(BUILD)/init.o: $(BUILD)/optField.o

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
$(BUILD)/output.o: $(BUILD)/sizeof.o
$(BUILD)/output.o: $(BUILD)/ice2_sub.o
$(BUILD)/output.o: $(BUILD)/optField.o

# finish.f90
$(BUILD)/finish.o: $(BUILD)/types.o

# ice.f90
$(BUILD)/ice.o: $(BUILD)/types.o
$(BUILD)/ice.o: $(BUILD)/atmosphere.o

# ice2.f90
$(BUILD)/ice2.o: $(BUILD)/types.o
$(BUILD)/ice2.o: $(BUILD)/atmosphere.o
$(BUILD)/ice2.o: $(BUILD)/update.o

# cleaning
TEMP = $(BUILD)/*.o $(BUILD)/*.mod $(BIN)/pinc
clean:
	rm -f $(TEMP)
