# 
# this Makefile compiles programs in "aux" directory
#

.SUFFIXES:

.PHONY: clean vclean

MODEL_E_ROOT = ..
MODEL_DIR = $(MODEL_E_ROOT)/model
DECKS_DIR = $(MODEL_E_ROOT)/decks
SCRIPTS_DIR = $(MODEL_E_ROOT)/exec

RUN     =
RUNDECK = $(DECKS_DIR)/$(RUN).R

OBJ_DIR = .
BIN_DIR = $(DECKS_DIR)/$(RUN)_bin

GCMLIB = $(BIN_DIR)/$(RUN).a

# getting options from configuration files
# default modelE configuration file:
include $(SCRIPTS_DIR)/modelErc
# configuration file with user-defined options:
sinclude ~/.modelErc         


# list of targets to be built by default
all: $(BIN_DIR)/qc $(BIN_DIR)/CMPE001 $(BIN_DIR)/conv_rsf $(BIN_DIR)/conv_gic $(BIN_DIR)/rsf2gic $(BIN_DIR)/conv_strat_rsf  $(BIN_DIR)/RMS $(BIN_DIR)/acc2nc
qflux: $(BIN_DIR)/vertflux.exe $(BIN_DIR)/ocnIC.exe $(BIN_DIR)/otspec.E001.exe 
# $(BIN_DIR)/rsfE001.addO.exe obsolete
odeep: $(BIN_DIR)/mkdeep.exe

# include file with implicit rules
include $(MODEL_DIR)/Rules.make

# look for .mod files in MODEL_DIR
ifeq ($(COMPILER),Absoft)
FFLAGS += -p $(MODEL_DIR)
else
FFLAGS += -I$(MODEL_DIR)
endif

#
#  targets
#

$(BIN_DIR)/CMPE001: CMPE001.f $(GCMLIB)
	$(F90) $(FFLAGS) $(CPPFLAGS) $(LFLAGS) -o $@ $^

$(BIN_DIR)/RMS: RMS.f $(GCMLIB)
	$(F90) $(FFLAGS) $(CPPFLAGS) $(LFLAGS) -o $@ $^

$(BIN_DIR)/conv_rsf: conv_rsf.o $(GCMLIB)
	$(F90) $(FFLAGS) $(CPPFLAGS) $(LFLAGS) -o $@ $^ $(LIBS)

conv_rsf.o: $(GCMLIB)

$(BIN_DIR)/conv_strat_rsf: conv_strat_rsf.o $(GCMLIB)
	$(F90) $(FFLAGS) $(CPPFLAGS) $(LFLAGS) -o $@ $^ $(LIBS)

conv_strat_rsf.o: $(GCMLIB)

$(BIN_DIR)/conv_gic: conv_gic.o $(GCMLIB)
	$(F90) $(FFLAGS) $(CPPFLAGS) $(LFLAGS) -o $@ $^ $(LIBS)

$(BIN_DIR)/rsf2gic: rsf2gic.o $(GCMLIB)
	$(F90) $(FFLAGS) $(CPPFLAGS) $(LFLAGS) -o $@ $^ $(LIBS)

conv_gic.o: $(GCMLIB)

rsf2gic.o: $(GCMLIB)

$(BIN_DIR)/qc: qc.o $(GCMLIB)
	$(F90) $(FFLAGS) $(CPPFLAGS) $(LFLAGS) -o $@ $^

qc.o: $(GCMLIB)

$(BIN_DIR)/vertflux.exe: vertflux.E001.o MAPS.o $(GCMLIB)
	$(F90) $(FFLAGS) $(CPPFLAGS) $(LFLAGS) -o $@ $^

vertflux.E001.o: $(GCMLIB)

$(BIN_DIR)/ocnIC.exe: ocnIC.E001.o MAPS.o $(GCMLIB)
	$(F90) $(FFLAGS) $(CPPFLAGS) $(LFLAGS) -o $@ $^

ocnIC.E001.o: $(GCMLIB)

$(BIN_DIR)/otspec.E001.exe: otspec.E001.o MAPS.o $(GCMLIB)
	$(F90) $(FFLAGS) $(CPPFLAGS) $(LFLAGS) -o $@ $^

otspec.E001.o: $(GCMLIB)

#$(BIN_DIR)/rsfE001.addO.exe: rsfE001.addO.o $(GCMLIB)
#	$(F90) $(FFLAGS) $(LFLAGS) -o $@ $^
#
#rsfE001.addO.o: $(GCMLIB)

$(BIN_DIR)/mkdeep.exe: mkdeep.o $(GCMLIB)
	$(F90) $(FFLAGS) $(CPPFLAGS) $(LFLAGS) -o $@ $^

mkdeep.o: $(GCMLIB)

$(GCMLIB): FORCE
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN)

ifdef NETCDFHOME
$(BIN_DIR)/acc2nc: acc2nc.o $(GCMLIB)
	$(F90) $(CPPFLAGS) $(LFLAGS) -o $@ $^ $(LIBS)
else
$(BIN_DIR)/acc2nc:
	@echo "NetCDF library is not specified"
	@echo "will skip $(BIN_DIR)/acc2nc"
endif

acc2nc.o: $(GCMLIB)

clean:
	rm -f *.LST *~ *.ERR *.CHK *.L

vclean:
	rm -f *.o *.mod *.vo *.inc .depend* *.sig

