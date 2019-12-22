# Makefile for Model E Atmosphere GCM B399+ 
.SUFFIXES:
.SUFFIXES: .o .U .S .GCM .OCN .f .F
.SUFFIXES:

.PHONY: rundeck depend gcm setup clean vclean newstart exe cmp .R htmldoc

#
# Print Usage if no targets were given
#
print_usage help all:
	@echo
	@echo '               GCM Model E Makefile'
	@echo 
	@echo 'This Makefile works with GNU make ("gmake" on most systems).'
	@echo 'It WILL NOT work with standard "make" which comes with the '
	@echo 'system (unless you are on Linux). If you see error messages'
	@echo 'above or below this message, that means that you are using'
	@echo 'wrong (not GNU) "make" command.'
	@echo 'This Makefile automatically recognizes and supports three'
	@echo 'types of systems: SGI, IBM workstation and Linux PC. Though'
	@echo '"setup" currently works only on SGIs.'
	@echo
	@echo The following is a list of supported targets:
	@echo
	@echo ' rundeck - create new rundeck'
	@echo ' depend  - create dependencies for specified rundeck'
	@echo ' gcm     - compile object files and build executable for specified rundeck'
	@echo ' aux     - compile standard auxiliary programs'
	@echo ' auxqflux- compile auxiliary programs for computing qflux'
	@echo ' auxdeep - compile auxiliary programs for setting deep ocean'
	@echo ' setup   - do setup for specified rundeck'
	@echo ' clean   - remove listings and error files'
	@echo ' vclean  - remove object files, .mod files and dependencies'
	@echo ' newstart- remove all files in the run directory'
	@echo ' exe     - compile gcm and put executable into RUN directory'
	@echo ' htmldoc - create web-based documentation for this RUN'
	@echo ' setup_nocomp - do not recompile, run only setup script'
	@echo
	@echo '          Recommended usage:'
	@echo 'The targets are made in order: "rundeck", "depend", "gcm", "setup",'
	@echo 'but since Makefile automatically creates necessary targets, you'
	@echo 'do not have to make all of them manually.'
	@echo 'You should always start with "rundeck", it will NOT be created automatically'
	@echo 'as it was in previous version. You may (and should) always skip "depend".'
	@echo 'So if you want to compile gcm and set up a run with a name "my_run"'
	@echo 'typically you will have to execute the commands:'
	@echo
	@echo '       gmake rundeck RUN=my_run'
	@echo '       gmake gcm RUN=my_run'
	@echo '       gmake setup RUN=my_run'
	@echo
	@echo 'You can actually skip "gmake gcm..." either and gcm will be'
	@echo 'compiled automatically, but you may not want to do it in most'
	@echo 'cases, since you want to look at the warning massages before'
	@echo 'doing "setup".'
	@echo
	@echo '          Customization'
	@echo 'You may choose not to see all the "technical" messages which'
	@echo 'gmake outputs to the screen. You can do it with option "-s"'
	@echo 'If you type:'
	@echo '        gmake -s setup RUN=my_run'
	@echo 'the ony information you will see will be compilation progress,'
	@echo 'warnings and errors.'
	@echo 'By default all warnings and errors are written to files with'
	@echo 'the name "fortran_file_name.ERR". If you prefer them to be'
	@echo 'sent to the screen you can do it by setting an environment'
	@echo 'variable MODELE_MAKE_OUTPUT=SCREEN .'
	@echo



#
# Directory structure and file names
#

MODEL_E_ROOT = ..
MODEL_DIR = .
DECKS_DIR = .
SCRIPTS_DIR = $(MODEL_E_ROOT)/exec

RUNSRC  = E001
RUN     =
RUNDECKSRC = $(MODEL_DIR)/$(RUNSRC).R
RUNDECK = $(DECKS_DIR)/$(RUN).R

OBJ_DIR = .
BIN_DIR = $(DECKS_DIR)/$(RUN)_bin
DOC_DIR = $(DECKS_DIR)/$(RUN)_doc

DEPENDFILE = $(OBJ_DIR)/.depend.$(RUN)

RUN_H = rundeck_opts.h
RUN_H_TMP = tmp_cpp.h

# getting options from configuration files
# default modelE configuration file:
include $(SCRIPTS_DIR)/modelErc
# configuration file with user-defined options:
sinclude ~/.modelErc

#
# Lists of objects and dependencies
# 

ifneq ($(MAKECMDGOALS),rundeck)
ifeq ($(shell [ -s $(RUNDECK) ] && echo yes),yes)

# extracting the list of objects (including options) from the rundeck
OBJ_LIST_O = $(shell [ -s $(RUNDECK) ] && perl \
	-e 'while(<> !~ /^Object +modules/ ){};' \
	-e 'while(($$_=<>) !~ /(^Data +input)|(^ *$$)/)' \
	-e '{chop; s/ *(!.*)*$$//; print " $$_";}' \
	< $(RUNDECK))

# stripping options from OBJ_LIST_O
OBJ_LIST = $(shell perl \
	-e '$$_="$(OBJ_LIST_O)"; s/\|[^|]*\|//g; print;' \
	)

# RFLAGS returns rundeck options for the current object (i.e. for $*)
RFLAGS = $(shell perl \
	-e '$$_="$(OBJ_LIST_O)"; m/\b$* *\|([^|]*)\|/; print " $$1";' \
	)

sinclude $(DEPENDFILE)
FSRCS = $(OBJ_LIST:=.f)
FSRCS_CPP = $(OBJ_LIST:=.f.cpp)
OBJS = $(OBJ_LIST:=.o)
# OBJSA - list of objects to include into gcmlib.a
OBJSA = $(OBJS)
ifneq ($(NETCDFHOME)blank,blank)
OBJSA += NCACC.o
endif
# extracting list of modules from from DEPENDFILE
MODS = $(shell perl \
	-e 'while(<>){ push @a, /(\S+\.mod)\s*:/gi; } print join " ",@a;' \
	$(DEPENDFILE) )
endif
endif


#
#include rules here
#

sinclude $(MODEL_DIR)/Rules.make


#########             TARGETS               ########

# Make a current rundeck
rundeck:
	@if [ ! -s $(RUNDECKSRC) ]; then \
	echo 'Source rundeck $(RUNDECKSRC) does not exist'; exit 1; fi
	@if [ "$(RUN)" = "" ]; then \
	echo 'You should specify run name on a command line (i.e. RUN=...)'; \
	exit 1; fi
	@echo '---- Regenerating Rundeck $(RUNDECK) from default ($(RUNDECKSRC)) ----'
	@if [ $(RUNDECK) != $(RUNDECKSRC) ]; then  \
	sed 's/^'$(RUNSRC)'\([:\. ]\)/'$(RUN)'\1/g' $(RUNDECKSRC) > $(RUNDECK); \
	fi;
	@echo '---------             Rundeck  created               ---------' 

#$(DEPENDFILE): depend

depend $(DEPENDFILE): $(RUNDECK) $(FSRCS) $(RUN_H)
ifeq ($(MACHINE),not_specified)
	@echo '*****      This architecture is not supported yet        *****'
	@echo '*****    Your architecture was detected as: $(UNAME)     *****'
	@echo '***** Corresponding section shoud be added to Rules.make *****'
	exit 1
endif
	@echo '--------          Rebuilding Dependencies            ---------'
	@[ -f $(RUN_H) ] || touch $(RUN_H)
	@echo "# dependencies generated by cpp:" > $(DEPENDFILE)
	@echo 'running CPP'
ifeq ($(MACHINE),IBM)
	for i in $(FSRCS); do $(CPP) -M $(CPPFLAGS) $$i > $$i.cpp; \
	cat $${i%.f}.u >> $(DEPENDFILE); done
else
ifeq ($(MACHINE),Linux)
# hack to get around use of 'traditional' in cpp call
	@for i in $(FSRCS); do $(CPP) -M $(CPPFLAGS) $$i | sed 's/\.f\.o/.o/g' >> $(DEPENDFILE); \
	$(CPP) $(CPPFLAGS) $$i > $$i.cpp; done
else
	@for i in $(FSRCS); do $(CPP) -M $(CPPFLAGS) $$i >> $(DEPENDFILE); \
	$(CPP) $(CPPFLAGS) $$i > $$i.cpp; done
endif
endif
	@echo 'running fmakedep'
	@$(FMAKEDEP) $(INCS) -f $(DEPENDFILE)_tmp $(FSRCS_CPP) 
	@cat $(DEPENDFILE)_tmp >> $(DEPENDFILE)
	@rm -f $(FSRCS_CPP) $(DEPENDFILE)_tmp
	@echo '---------         done with dependencies             ---------'

$(RUN_H): $(RUNDECK) FORCE
	-@rm -f $(RUN_H_TMP)
	@perl -e 'while(($$_=<>) && ! /^ *Preprocessor +Options/i ){};' \
	-e 'while(($$_=<>) && ! /^ *End +Preprocessor +Options/i )' \
	-e '{chop; s/ *(!.*)*$$//; $$_ && print "$$_\n";}' \
	< $(RUNDECK) > $(RUN_H_TMP)
	@[ -f $(RUN_H) ] && diff $(RUN_H_TMP) $(RUN_H) >/dev/null 2>&1  \
	|| ( mv $(RUN_H_TMP) $(RUN_H); echo 'Updated CPP Options' )

$(BIN_DIR):
	mkdir $(BIN_DIR)

# Create an executable
gcm $(BIN_DIR)/$(RUN).exe: $(RUNDECK) $(OBJS)
	@echo linking executable $(MSG)
	$(F90) $(LFLAGS) $(OBJS) $(LIBS) -o $(RUN).exe  $(LINK_OUTPUT)
ifdef LINK_OUTPUT
	@if [ -s $(RUN).ERR ] ; then cat $(RUN).ERR; else rm -f $(RUN).ERR; fi
endif
	@if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	mv $(RUN).exe $(BIN_DIR)
	@echo '--------------------------------------------------------------'
	@echo '---------       GCM successfuly compiled             ---------'
	@echo '---------    executable $(RUN).exe was created     ---------'
	@echo '--------------------------------------------------------------'

gcmlib $(BIN_DIR)/$(RUN).a: $(RUNDECK) $(MODS) $(OBJSA)
	ar cr $(RUN).a $(OBJSA)
	@if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	mv $(RUN).a $(BIN_DIR)


# Setup the gcm
setup:  $(BIN_DIR)/$(RUN).exe
	@echo '---------        Starting setup for $(RUN)         ---------'
	@echo '--------------------------------------------------------------'
	@$(SETUP) $(RUN)

# Mkexe the gcm
exe:	$(BIN_DIR)/$(RUN).exe
	@if [ ! -s /u/cmrun/$(RUN)/I ] ; then \
	echo 'directory /u/cmrun/$(RUN) or input files inside it are missing';\
	echo 'you should make "setup" first'; \
	exit 1; fi
	cp -f $(BIN_DIR)/$(RUN).exe /u/cmrun/$(RUN)/

newstart:
	@if [ "$(RUN)" = "" ]; then \
        echo 'You should specify run name on a command line (i.e. RUN=...)'; \
        exit 1; fi
	rm -i /u/cmrun/$(RUN)/*

htmldoc: $(FSRCS_CPP)
	[ -d $(DOC_DIR) ] || mkdir $(DOC_DIR)
	$(SCRIPTS_DIR)/gcmdoc.pl -O $(DOC_DIR) -R $(RUN) -CPP f.cpp $(FSRCS)
	rm -f $(FSRCS_CPP)

clean:
	rm -f *.LST *~ *.ERR *.CHK *.L *.cpp

vclean:
	rm -f *.o *.mod *.vo *.inc .depend* *.sig

src_cpp_dir: $(FSRCS_CPP)
	mkdir $(DECKS_DIR)/$(RUN)_src
	mv $(FSRCS_CPP)  $(DECKS_DIR)/$(RUN)_src
	cd $(DECKS_DIR)/$(RUN)_src; ls `pwd`/*.cpp >file_list

