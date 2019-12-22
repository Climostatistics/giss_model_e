#
# this file contains rules shared by Makefiles in "model" and "aux"
#

.PHONY: FORCE


######  Some user customizable settings:   ########

# EXTRA_FFLAGS specifies some extra flags you want to pass 
# to Fortarn compiler, like
# -g        - include debugging information
# -listing  - create listings (.L)
EXTRA_FFLAGS =

# CPPFLAGS specifies some pre-processor directives you want to pass 
# to Fortarn compiler, like
# -DCOMPILER=Absoft    - compile code that is specific to Absoft
CPPFLAGS =

ifeq ($(OUTPUT_TO_FILES),YES)
  COMP_OUTPUT = > $*.ERR 2>&1 || { r=$$? ; cat $*.ERR ; exit $$r ; }
  LINK_OUTPUT = > $(RUN).ERR 2>&1 || { r=$$? ; cat $(RUN).ERR ; exit $$r ; }
else
  COMP_OUTPUT =
  LINK_OUTPUT =
endif

# if -s specified enable some extra messages
ifeq ($(findstring s,$(MFLAGS)),s)
MSG = 
else
MSG = > /dev/null
endif

#
# starting machine - specific options
#

NO_COMMAND = echo Requested target is not supported on $(UNAME); exit 1;
F90 = $(NO_COMMAND)
FMAKEDEP = $(NO_COMMAND)
F =  $(NO_COMMAND)
U = $(NO_COMMAND)
SETUP = $(SCRIPTS_DIR)/setup_e.pl
SETUP_GFDL = $(SCRIPTS_DIR)/setup_e_gfdl.pl
CPP = $(NO_COMMAND)
MACHINE = not_specified
LIBS =
INCS =
F90_VERSION = 'Unknown compiler version'

UNAME = $(shell uname)


# SGI - specific options here
ifeq ($(UNAME),IRIX64)
MACHINE = SGI
F90 = f90
CPP = /lib/cpp -P
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend -H
F       = $(SCRIPTS_DIR)/fco2_90
U	= $(SCRIPTS_DIR)/uco2_f90
CPPFLAGS = -DMACHINE_SGI
#FFLAGS = -cpp -O2 -64 -mips4 -OPT:reorg_comm=off -w2 -OPT:Olimit=5745
FFLAGS = -ftpp -O2 -64 -mips4 -OPT:reorg_comm=off -w2 -OPT:Olimit=6000
LFLAGS = -64 -O2 -mips4 -lfastm -OPT:reorg_common=OFF
ifeq ($(MP),YES)
FFLAGS += -mp
LFLAGS += -mp
endif
# suppress some linker warnings if no verbose output
ifeq ($(VERBOSE_OUTPUT),NO)
LFLAGS += -LD_MSG:OFF=84,85,15,134
endif
# uncomment next two lines for extra debugging
#FFLAGS += -DEBUG:div_check=3 -DEBUG:subscript_check=ON -DEBUG:trap_uninitialized=ON
#LFLAGS += -DEBUG:conform_check=YES -DEBUG:div_check=3 -DEBUG:subscript_check=ON -DEBUG:trap_uninitialized=ON
# not sure if the following will help the debugging ...
# FFLAGS += -DEBUG:verbose_runtime=ON
# LFLAGS += -DEBUG:verbose_runtime=ON
F90_VERSION = $(shell $(F90) -version 2>&1)
endif

# SGI-32 - specific options here
ifeq ($(UNAME),IRIX)
MACHINE = SGI
F90 = f90
CPP = /lib/cpp -P
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend -H
F       = $(SCRIPTS_DIR)/fco2_90
U       = $(SCRIPTS_DIR)/uco2_f90
CPPFLAGS = -DMACHINE_SGI
FFLAGS = -cpp -O2 -mips4 -OPT:reorg_comm=off -w2 -OPT:Olimit=5745
LFLAGS = -O2 -mips4 -lfastm -mp -OPT:reorg_common=OFF -Wl,-woff,134 -Wl,-woff,15
F90_VERSION = $(shell $(F90) -version 2>&1)
endif


# Linux - specific options here
# If COMPILER is set will use options for that compiler, otherwise will
# use the default (Absoft). The following compilers are recognized:
# Absoft, Lahey, PGI, Intel (not working yet), Vast (not working yet).
ifeq ($(UNAME),Linux)
MACHINE = Linux
ifndef COMPILER
COMPILER = Absoft
endif

## This is for the VAST compiler (DOES NOT WORK)
ifeq ($(COMPILER),Vast)
F90 = f90
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend -m vo
endif

## this is for Intel Compiler (DOES NOT WORK WELL)
ifeq ($(COMPILER),Intel)
F90 = ifc
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend -m d -int
FFLAGS = -O2
LFLAGS = -O2
endif

## This is for the Lahey/Fujitsu compiler
ifeq ($(COMPILER),Lahey)
F90 = lf95
CPP = /usr/bin/cpp -P -traditional
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend 
CPPFLAGS = -DCONVERT_BIGENDIAN -DMACHINE_Linux
FFLAGS = -O
LFLAGS = 
F90_VERSION = $(shell $(F90) --version)
endif

## This is for the Absoft PROfortran compiler
ifeq ($(COMPILER),Absoft)
F90 = f90
CPP = /usr/bin/cpp -P -traditional
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend -h
CPPFLAGS = -DCONVERT_BIGENDIAN -DMACHINE_Linux
FFLAGS = -O2
LFLAGS = -lf90math -lV77 -lU77
# uncomment next two lines for extensive debugging
#FFLAGS += -g -trap=INVALID,DIVBYZERO,OVERFLOW -B111
#LFLAGS += -lefence -g
endif

## This is for the Portland Group Compiler (PGI)
ifeq ($(COMPILER),PGI)
F90 = pgf90 -Mbyteswapio
CPP = /usr/bin/cpp -P -traditional
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
CPPFLAGS = -DMACHINE_Linux
FFLAGS = -O2
LFLAGS = 
ifeq ($(MP),YES)   # edit for PGI OpenMP compatability???
FFLAGS += -Msmp
LFLAGS += -Msmp
endif
# uncomment next two lines for extensive debugging
#FFLAGS += 
#LFLAGS += 
endif
endif

# IBM - specific options here
ifeq ($(UNAME),AIX)
MACHINE = IBM
F90 = xlf90_r
CPP = /lib/cpp -P
FMAKEDEP = perl $(SCRIPTS_DIR)/sfmakedepend
# ibm compiler doesn't understand "-D" . Have to use "-WF,-D..."
CPPFLAGS =
FFLAGS = -O2 -qfixed -qsuffix=cpp=f -qmaxmem=16384 -WF,-DMACHINE_IBM
# one may need to add -bmaxstack:0x1000000 if rusns out of stack
LFLAGS = -O2 -bmaxdata:0x10000000
# no guarantee that the following line gives correct info
F90_VERSION = $(shell what /usr/lpp/xlf/bin/xlfentry | tail -1)
endif


# DEC Alpha - specific options here
# This option has not been tested yet. It may need some adjustment.
ifeq ($(UNAME),OSF1)
MACHINE = DEC
F90 = f90
CPP = /lib/cpp -P
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
# CPPFLAGS = -DCONVERT_BIGENDIAN -DMACHINE_DEC
# FFLAGS = -O2 -cpp
# LFLAGS = -O2
CPPFLAGS = -DMACHINE_DEC
FFLAGS = -O2 -cpp -convert big_endian
LFLAGS = -O2 -convert big_endian
ifeq ($(MP),YES)
FFLAGS += -omp
LFLAGS += -omp
endif
F90_VERSION = $(shell $(F90) -version 2>&1)
endif


# end of machine - specific options

#
# Check for extra options specified in modelErc
#

ifdef NETCDFHOME
ifeq ($(MACHINE),SGI)
  LIBS += -L$(NETCDFHOME)/lib64 -lnetcdf
else
  LIBS += -L$(NETCDFHOME)/lib -lnetcdf
endif
  FFLAGS += -I$(NETCDFHOME)/include
  INCS += -I $(NETCDFHOME)/include
endif

#
# Pattern  rules
#

FORCE:

*.mod: FORCE

%.mod: 
	@echo -n checking $@:
#	@echo 'called rule for $@, depends on $^ built from: '`cat $@.sig`
	@if [ "$<empty" = "empty" ]; then \
	echo "No dependency for $@ : ";\
	echo "Unsupported architecture or Makefile error";\
	exit 1; fi
	@if [ ! -s $@.sig ] || [ "`cat $@.sig`" != "$<" ]; then \
	echo "  dependencies for $@ have changed - recompiling"; \
	rm -f $<; $(MAKE) $< RUN=$(RUN); else echo '  ok'; fi

# Standard fortran
# .timestemp is a hack to set proper times on .o and .mod
# For the Absoft/Lahey/PGI compilers, we need to force a cpp run through
%.o: %.f
	@echo -n compiling $< ... $(MSG)
	@touch .timestamp
ifeq ($(COMPILER),Absoft)
	cp $*.f $*.F
	$(F90) -c $(FFLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS) $(RFLAGS) $*.F \
	  $(COMP_OUTPUT)
	rm -f $*.F
else
ifeq ($(COMPILER),Lahey)
	cp $*.f $*.F
	$(F90) -c $(FFLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS) $(RFLAGS) $*.F \
	  $(COMP_OUTPUT)
	rm -f $*.F
else
ifeq ($(COMPILER),PGI)
	cp $*.f $*.F
	$(F90) -c $(FFLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS) $(RFLAGS) $*.F \
	  $(COMP_OUTPUT)
	rm -f $*.F
else
	$(F90) -c $(FFLAGS) $(EXTRA_FFLAGS) $(CPPFLAGS) $(RFLAGS) $*.f \
	  $(COMP_OUTPUT)
endif
endif
endif
	-@if [ `ls | grep ".mod" | tail -1` ] ; then for i in *.mod; \
	  do if [ ! -s $$i.sig ] || [ `find $$i -newer $$i.sig` ] ; then \
	  echo $@ > $$i.sig; fi; done; fi 
	@touch -r .timestamp $@
	@if [ -s $*.ERR ] ; then echo $(MSG); else echo Done $(MSG); fi
ifdef COMP_OUTPUT
	@if [ -s $*.ERR ] ; then cat $*.ERR; else  rm -f $*.ERR; fi
endif

# cpp preprocessing
%.f.cpp: %.f
	@echo preprocessing $<  $(MSG)
	$(CPP) $(CPPFLAGS) $*.f > $*.f.cpp

# Update files
%.o: %.U
	$(U)  $*

# GISS-Fortran source files (with line numbers) 
%.o: %.S
	$(F) $<

# MAPS.GCM
%.o: %.GCM
	$(F) $<

# FUNTABLE.OCN
%.o: %.OCN
	$(F) $<


# end of Pattern  rules

