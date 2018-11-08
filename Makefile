# Makefile of pinc_MPI_bg

FC = scorep --nomemory mpif90
FCFLAGS=-O3 -real-size 64 -traceback -unroll=4 -ip
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
	pinc.o

# add build directory as prefix to path of %.o files
OBJ=$(addprefix $(BUILD)/, $(OFILES))

# general rules
$(BUILD)/%.o: $(SOURCE)/%.f90
	$(FC) $(FCFLAGS) -module $(BUILD) -c $< -o $@

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



# program for file difference
DIFFOBJ = types.o \
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
l2diff:$(DIFFOBJ)
	$(FC) $(FCFLAGS) -o l2diff $(DIFFOBJ)

# programme for making tec360 layout files
LAYOBJ = makeLayout.o types.o
layout: $(LAYOBJ)
	$(FC) $(FCFLAGS) -o makeLayout $(LAYOBJ)

# test algebra_module
algebra: algebra.mod
	$(FC) $(FCFLAGS) -o testAlgebra testAlgebra.f90 algebra.f90

# test xweno_module
XOBJ = 	types.o xweno.o testXWENO.o debug.o
xweno:	$(XOBJ)
	$(FC) $(FCFLAGS) -o testXWENO $(XOBJ)

# cleaning
TEMP = $(BUILD)/*.o $(BUILD)/*.mod $(BIN)/pinc
clean:
	rm -f $(TEMP)

## Loesche Tecplot files
#cleandat:
#	rm -f *.dat *.dat~ *.res
#
#
## transform all DAT-files to PLT-Files
#
#
#DATFILES = $(shell ls *.dat)
#PLTFILES = $(DATFILES:%.dat=%.plt)
#
## Hier schien mir der Doppelpunkt bei der oberen Zeile zu viel. 
## Wenn der Eintrag mit der shell in der Klammer tut (hab ich noch nie
## probiert), muesste es klappen
#
#
#plt: $(PLTFILES)
#
## Hier hab ich die leere Anweisung wegoperiert. 
## Die Abhaengigkeit muesste reichen. 
#
## Pattern rule fuer die Umwandlung dat -> plt
#%.plt : %.dat
#	preplot $*.dat
#
#
## mit dem preplot Aufruf bin ich mir nicht ganz sicher, ob der automatisch 
## die Endung .plt macht. Es ist einfach zu lange her, dass ich tecplot 
## zur Verfuegung hatte.
