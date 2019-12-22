# Makefile for Model II' diagnostic output routines
# Begin macro definitions

F	= f90
FOPTS	= -64 -O2
NCLIB	= /usr/local/netcdf-3.4/lib64/libnetcdf.a

IJDAG	= ijdag
JLDAG	= jldag
IJKDAG	= ijkdag
LATDAG	= latdag

PROGS	= $(IJDAG) $(JLDAG) $(IJKDAG) $(LATDAG)

NCINP	= ncinp.o
NCOUT	= ncout.o

UTLOBJ	= $(NCINP) $(NCOUT)

# Make rules
# first line is for turning off existing suffix rules
.SUFFIXES:
.SUFFIXES: .o .f

# Create object files from Fortran source files
.f.o:
	$(F) -c $(FOPTS) $<

# Create executable files from Fortran source files
.f:
	$(F) $(FOPTS) $< $(UTLOBJ) $(NCLIB) -o $*

# Make executable files dependent on io modules and source files

$(PROGS):	$(UTLOBJ)

# compile io utilities
utl:	$(UTLOBJ)

# compile all executables
all:	$(PROGS)

clean:
	rm -f *.o *.mod $(PROGS)
