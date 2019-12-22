.SUFFIXES:
.PHONY: rundeck depend gcm setup clean vclean newstart exe cmp .R htmldoc setup_script setup_nocomp

# modelE directory structure
MODEL_E_ROOT = ..
MODEL_DIR = $(MODEL_E_ROOT)/model
AUX_DIR = $(MODEL_E_ROOT)/aux
SCRIPTS_DIR = $(MODEL_E_ROOT)/exec
DECKS_DIR = $(shell pwd)
BIN_DIR = $(DECKS_DIR)/$(RUN)_bin

USER = $(shell whoami)
HOMEDIR = $(wildcard ~)
MODELERC = $(HOMEDIR)/.modelErc
# getting options from configuration files
# default modelE configuration file:
include $(SCRIPTS_DIR)/modelErc
# configuration file with user-defined options:
sinclude $(MODELERC)

ifeq ($(VERBOSE_OUTPUT),NO)
  MAKEFLAGS=-s
endif

# SETUP = $(SCRIPTS_DIR)/setup_e

RUNSRC = E001

ALIAS_LIST = CMPE001 qc

help:
	$(MAKE) -C $(MODEL_DIR) help

sinclude $(MODEL_DIR)/Rules.make

rundeck:
	@if [ -f $(DECKS_REPOSITORY)/$(RUN).R ] && \
	! find $(DECKS_REPOSITORY)/$(RUN).R -user $(USER) > /dev/null; then \
	echo 'name $(RUN) is already used by somebody else'; \
	echo 'please choose a unique one'; exit 1; fi
	@if [ -s $(DECKS_REPOSITORY)/$(RUN).R -a $(OVERWRITE) != YES ];\
	then echo '$(RUN).R already exists in repository'; \
	echo 'use OVERWRITE=YES to overwrite'; \
	exit 1; fi
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR) RUNSRC=$(RUNSRC)
	@cp -f $(RUN).R $(DECKS_REPOSITORY)/$(RUN).R
	@chmod 644 $(DECKS_REPOSITORY)/$(RUN).R      #just in case ...

gcm $(BIN_DIR)/$(RUN).exe depend:
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR)

setup: gcm
	$(MAKE) setup_script RUN=$(RUN)

setup_nocomp:
	@[ -d $(BIN_DIR) ] || mkdir $(BIN_DIR)
	@cp $(CMRUNDIR)/$(RUN)/$(RUN).exe $(BIN_DIR)/$(RUN).exe || \
	{ echo "No executable in run directory"; \
	echo "You have to run 'gmake setup ...' first"; exit 1; }
	$(MAKE) setup_script RUN=$(RUN)

setup_script:
	@echo '---------       Looks like it was compiled OK        ---------'
	@echo '----- Saving Rundeck and other info to global repository -----'
	@if [ ! -d $(DECKS_REPOSITORY) ]; then \
	  echo '--- $(DECKS_REPOSITORY) does not exist. Rundeck not saved ---';\
	else \
	  if [ -f $(DECKS_REPOSITORY)/$(RUN).R ] && \
	  ! find $(DECKS_REPOSITORY)/$(RUN).R -user $(USER) > /dev/null; then \
	    echo 'name $(RUN) is already used by somebody else'; \
	    echo 'please choose a unique one'; exit 1; \
	  else \
	    cp -f $(RUN).R $(DECKS_REPOSITORY)/$(RUN).R ;\
	    chmod 644 $(DECKS_REPOSITORY)/$(RUN).R  ;    \
	    echo '### Information below describes your run. Do not delete! ###'\
	    >> $(DECKS_REPOSITORY)/$(RUN).R ; \
	    date >> $(DECKS_REPOSITORY)/$(RUN).R ; \
	    echo $(F90_VERSION) >> $(DECKS_REPOSITORY)/$(RUN).R ;\
	    if [ -f $(MODEL_DIR)/CVS/Tag ]; then \
	      echo "CVS Repository Tag: `cat $(MODEL_DIR)/CVS/Tag`" \
	      >> $(DECKS_REPOSITORY)/$(RUN).R; else \
	      echo 'CVS Repository: MAIN Branch' \
	      >> $(DECKS_REPOSITORY)/$(RUN).R ;\
	    fi ;\
	  fi ;\
	fi
	@echo '---------        Starting setup for $(RUN)          ---------'
	@echo '--------------------------------------------------------------'
	@$(SETUP) $(RUN)

$(SETUP_GFDL):
	$(SCRIPTS_DIR)/create_setup_gfdl

setup_gfdl: gcm $(SETUP_GFDL)
	@echo '---------       Looks like it was compiled OK        ---------'
	@echo '----- Saving Rundeck and other info to global repository -----'
	@if [ ! -d $(DECKS_REPOSITORY) ]; then \
	  echo '--- $(DECKS_REPOSITORY) does not exist. Rundeck not saved ---';\
	else \
	  if [ -f $(DECKS_REPOSITORY)/$(RUN).R ] && \
	  ! find $(DECKS_REPOSITORY)/$(RUN).R -user $(USER) > /dev/null; then \
	    echo 'name $(RUN) is already used by somebody else'; \
	    echo 'please choose a unique one'; exit 1; \
	  else \
	    cp -f $(RUN).R $(DECKS_REPOSITORY)/$(RUN).R ;\
	    chmod 644 $(DECKS_REPOSITORY)/$(RUN).R  ;    \
	    echo '### Information below describes your run. Do not delete! ###'\
	    >> $(DECKS_REPOSITORY)/$(RUN).R ; \
	    date >> $(DECKS_REPOSITORY)/$(RUN).R ; \
	    echo $(F90_VERSION) >> $(DECKS_REPOSITORY)/$(RUN).R ;\
	    if [ -f $(MODEL_DIR)/CVS/Tag ]; then \
	      echo "CVS Repository Tag: `cat $(MODEL_DIR)/CVS/Tag`" \
	      >> $(DECKS_REPOSITORY)/$(RUN).R; else \
	      echo 'CVS Repository: MAIN Branch' \
	      >> $(DECKS_REPOSITORY)/$(RUN).R ;\
	    fi ;\
	  fi ;\
	fi
	@echo '---------        Starting setup for $(RUN)          ---------'
	@echo '--------------------------------------------------------------'
	@$(SETUP_GFDL) $(RUN)


aux:
	$(MAKE) -C $(AUX_DIR) RUN=$(RUN) DECKS_DIR=$(DECKS_DIR)

auxqflux: aux
	$(MAKE) qflux -C $(AUX_DIR) RUN=$(RUN) DECKS_DIR=$(DECKS_DIR)

auxdeep: aux
	$(MAKE) odeep -C $(AUX_DIR) RUN=$(RUN) DECKS_DIR=$(DECKS_DIR)

clean vclean:
	$(MAKE) -C $(MODEL_DIR) $@
	$(MAKE) -C $(AUX_DIR) $@

exe:  $(BIN_DIR)/$(RUN).exe
	@if ! find $(DECKS_REPOSITORY)/$(RUN).R -user $(USER) >/dev/null; then\
	echo 'Either you didn''t run setup yet or the rundeck name $(RUN)'; \
	echo 'belongs to somebody else'; echo 'please check'; exit 1; fi
	@if [ ! -s $(CMRUNDIR)/$(RUN)/I ] ; then \
	echo \
	  'directory $(CMRUNDIR)/$(RUN) or input files inside it are missing';\
	echo 'you should make "setup" first'; \
	exit 1; fi
	mv -f $(BIN_DIR)/$(RUN).exe $(CMRUNDIR)/$(RUN)/
	@echo 'adding corresponding info to rundeck repository'
	@perl -e "while(<>){ if(/^### Information/){exit;} print; }" \
	<$(DECKS_REPOSITORY)/$(RUN).R  >$(DECKS_REPOSITORY)/$(RUN).R.tmp
	@echo '### Executable was recompiled ' >> $(DECKS_REPOSITORY)/$(RUN).R
	@date >> $(DECKS_REPOSITORY)/$(RUN).R
	@echo $(F90_VERSION) >> $(DECKS_REPOSITORY)/$(RUN).R
	@if [ -f $(MODEL_DIR)/CVS/Tag ]; then \
	echo "CVS Repository Tag: `cat $(MODEL_DIR)/CVS/Tag`" \
	>> $(DECKS_REPOSITORY)/$(RUN).R; else \
	echo 'CVS Repository: MAIN Branch' \
	>> $(DECKS_REPOSITORY)/$(RUN).R; fi
	@echo 'diff old.R new.R :' >> $(DECKS_REPOSITORY)/$(RUN).R
	@if ! diff $(DECKS_REPOSITORY)/$(RUN).R.tmp $(RUN).R \
	>> $(DECKS_REPOSITORY)/$(RUN).R; then echo 'rundeck changed'; fi
	@rm -f $(DECKS_REPOSITORY)/$(RUN).R.tmp
	@echo done

newstart:
	@if [ "$(RUN)" = "" ]; then \
	echo 'You should specify run name on a command line (i.e. RUN=...)'; \
	exit 1; fi
	rm -i $(CMRUNDIR)/$(RUN)/*

alias: aux
	for i in $(ALIAS_LIST); do ln -sf $(BIN_DIR)/$$i .; done

unalias:
	for i in $(ALIAS_LIST); do rm -f $$i; done

htmldoc:
	$(MAKE) -C $(MODEL_DIR) $@ RUN=$(RUN) DECKS_DIR=$(DECKS_DIR)

config configure:
	@if [ -f $(MODELERC) -a $(OVERWRITE) != YES ]; then \
	echo '~/.modelErc already exists, use OVERWRITE=YES to overwrite'; \
	exit 1; fi
	@cp -f $(SCRIPTS_DIR)/modelErc $(MODELERC)
	@if [ `hostname` = *'giss.nasa.gov' ]; then \
	sed "s/MAILTO=/MAILTO=`whoami`@giss.nasa.gov/" < $(MODELERC) > \
	$(MODELERC).tmp; mv -f $(MODELERC).tmp $(MODELERC); fi
	@echo '$(MODELERC) created.'
	@echo 'You may want to edit it manually to customize the options.'

update:
	$(MAKE) -C $(MODEL_E_ROOT) $@ RELEASE=$(RELEASE)


