# Makefile of pinc_MPI_bg

FC = mpif90
FCFLAGS=-O3 -real-size 64 -traceback -unroll=4 -ip
#FCFLAGS=-O0 -fpe0 -check bounds -real-size 64 -traceback -unroll=4 -ip
LIPNAG =
LIBHYPRE = /home/atmodynamics/achatz/Installation/hypre-2.11.2/src/lib

OBJ =	types.o \
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

# general rules
%.o: %.f90 
	$(FC) $(FCFLAGS) -c $<

# the main target
pinc:$(OBJ)
	$(FC) $(FCFLAGS) -o pinc $(OBJ) $(LIPNAG77) -L$(LIBHYPRE) -lHYPRE
#-------------------
#  dependencies
#-------------------

# fluxes.f90
fluxes.o: types.o
fluxes.o: xweno.o
fluxes.o: muscl.o
fluxes.o: atmosphere.o
fluxes.o: algebra.o


# xweno.f90
xweno.o: types.o
xweno.o: debug.o

# debug.f90
debug.o: types.o
debug.o: atmosphere.o

# poisson.f90
poisson.o: types.o
poisson.o: mpi.o

mpi.o: fluxes.o
mpi.o: types.o
timeScheme.o: types.o
atmosphere.o: types.o
init.o: types.o
muscl.o: types.o
wkb.o: types.o
boundary.o: types.o
update.o: types.o
output.o: types.o
finish.o: types.o



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
TEMP = *.o *.mod pinc
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
