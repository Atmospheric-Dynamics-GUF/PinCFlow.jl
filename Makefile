# Makefile of PinCFlow

# Set the compiler.
FC = mpif90
# FC = mpiifort

# Set the debug option.
DEBUGFLAG = false
# DEBUGFLAG = true

# Set flags.
ifeq ($(FC), mpiifort) # Intel
  ifeq ($(DEBUGFLAG), true)
    FCFLAGS = -O0 -fpe0 -check all -real-size 64 -traceback -unroll=4 -ip -debug full
  else ifeq ($(DEBUGFLAG), false)
    FCFLAGS = -O3 -real-size 64 -traceback -unroll=4 -ip
  endif
  MODULEFLAG = -module $(BUILD)
else ifeq ($(FC), mpif90) # GNU
  ifeq ($(DEBUGFLAG), true)
    FCFLAGS = -Og -g -fallow-argument-mismatch -fcheck=all -Wall -Wno-unused-variable -fdefault-real-8 -fbacktrace -funroll-loops -Wno-unused-dummy-argument -Wno-conversion-extra
  else ifeq ($(DEBUGFLAG), false)
    FCFLAGS = -O3 -fdefault-real-8 -fbacktrace -funroll-loops -fallow-argument-mismatch -w
  endif
  MODULEFLAG = -J$(BUILD)
endif
NCFLAGS = `nf-config --fflags --flibs` -Wl,-rpath,`nf-config --prefix`/lib

# Define directories for sources and binaries.
BIN = $(shell mkdir -p ./bin) ./bin
BUILD = $(shell mkdir -p ./build) ./build
SOURCE = ./src

# Define object files.
OFILES = atmosphere.o \
		bicgstab_tools.o \
		boundary.o \
		finish.o \
		fluxes.o \
		init.o \
		mpi.o \
		muscl.o \
		output_netcdf.o \
		pinc.o \
		poisson.o \
		sizeof.o \
		timeScheme.o \
		types.o \
		update.o

# Add build directory as prefix to path of object files.
OBJ = $(addprefix $(BUILD)/, $(OFILES))

# Set the rule for object files.
$(BUILD)/%.o: $(SOURCE)/%.f90
	$(FC) $(FCFLAGS) $(MODULEFLAG) -c $< -o $@ $(NCFLAGS)

# Set the rule for the executable.
pinc:$(OBJ)
	$(FC) $(FCFLAGS) -o $(BIN)/pinc $(OBJ) $(NCFLAGS)

# Set the rule for cleaning.
TEMP = $(BUILD)/*.o $(BUILD)/*.mod $(BIN)/pinc
clean:
	rm -f $(TEMP)

# List dependencies of atmosphere.f90.
$(BUILD)/atmosphere.o: $(BUILD)/sizeof.o
$(BUILD)/atmosphere.o: $(BUILD)/types.o

# List dependencies of bicgstab_tools.f90.
$(BUILD)/bicgstab_tools.o: $(BUILD)/atmosphere.o
$(BUILD)/bicgstab_tools.o: $(BUILD)/mpi.o
$(BUILD)/bicgstab_tools.o: $(BUILD)/timeScheme.o
$(BUILD)/bicgstab_tools.o: $(BUILD)/types.o

# List dependencies of boundary.f90.
$(BUILD)/boundary.o: $(BUILD)/atmosphere.o
$(BUILD)/boundary.o: $(BUILD)/fluxes.o
$(BUILD)/boundary.o: $(BUILD)/mpi.o
$(BUILD)/boundary.o: $(BUILD)/types.o

# List dependencies of finish.f90.
$(BUILD)/finish.o: $(BUILD)/types.o

# List dependencies of fluxes.f90.
$(BUILD)/fluxes.o: $(BUILD)/atmosphere.o
$(BUILD)/fluxes.o: $(BUILD)/muscl.o
$(BUILD)/fluxes.o: $(BUILD)/sizeof.o
$(BUILD)/fluxes.o: $(BUILD)/types.o

# List dependencies of init.f90.
$(BUILD)/init.o: $(BUILD)/atmosphere.o
$(BUILD)/init.o: $(BUILD)/boundary.o
$(BUILD)/init.o: $(BUILD)/mpi.o
$(BUILD)/init.o: $(BUILD)/output_netcdf.o
$(BUILD)/init.o: $(BUILD)/sizeof.o
$(BUILD)/init.o: $(BUILD)/types.o

# List dependencies of mpi.f90.
$(BUILD)/mpi.o: $(BUILD)/fluxes.o
$(BUILD)/mpi.o: $(BUILD)/types.o

# List dependencies of muscl.f90.
$(BUILD)/muscl.o: $(BUILD)/types.o

# List dependencies of output_netcdf.f90.
$(BUILD)/output_netcdf.o: $(BUILD)/atmosphere.o
$(BUILD)/output_netcdf.o: $(BUILD)/sizeof.o
$(BUILD)/output_netcdf.o: $(BUILD)/types.o

# List dependencies of pinc.f90.
$(BUILD)/pinc.o: $(BUILD)/atmosphere.o
$(BUILD)/pinc.o: $(BUILD)/bicgstab_tools.o
$(BUILD)/pinc.o: $(BUILD)/boundary.o
$(BUILD)/pinc.o: $(BUILD)/finish.o
$(BUILD)/pinc.o: $(BUILD)/fluxes.o
$(BUILD)/pinc.o: $(BUILD)/init.o
$(BUILD)/pinc.o: $(BUILD)/mpi.o
$(BUILD)/pinc.o: $(BUILD)/output_netcdf.o
$(BUILD)/pinc.o: $(BUILD)/poisson.o
$(BUILD)/pinc.o: $(BUILD)/sizeof.o
$(BUILD)/pinc.o: $(BUILD)/timeScheme.o
$(BUILD)/pinc.o: $(BUILD)/types.o
$(BUILD)/pinc.o: $(BUILD)/update.o

# List dependencies of poisson.f90.
$(BUILD)/poisson.o: $(BUILD)/atmosphere.o
$(BUILD)/poisson.o: $(BUILD)/bicgstab_tools.o
$(BUILD)/poisson.o: $(BUILD)/mpi.o
$(BUILD)/poisson.o: $(BUILD)/sizeof.o
$(BUILD)/poisson.o: $(BUILD)/timeScheme.o
$(BUILD)/poisson.o: $(BUILD)/types.o

# List dependencies of timeScheme.f90.
$(BUILD)/timeScheme.o: $(BUILD)/types.o

# List dependencies of tracer.f90.
$(BUILD)/tracer.o: $(BUILD)/atmosphere.o
$(BUILD)/tracer.o: $(BUILD)/types.o

# List dependencies of update.f90.
$(BUILD)/update.o: $(BUILD)/atmosphere.o
$(BUILD)/update.o: $(BUILD)/boundary.o
$(BUILD)/update.o: $(BUILD)/fluxes.o
$(BUILD)/update.o: $(BUILD)/mpi.o
$(BUILD)/update.o: $(BUILD)/poisson.o
$(BUILD)/update.o: $(BUILD)/timeScheme.o
$(BUILD)/update.o: $(BUILD)/types.o