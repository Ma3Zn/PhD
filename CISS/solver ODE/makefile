EXE_DIR := ./

LINKFLAGS 	:= -o
LINK_LIB	:= -lgsl -lgslcblas -lm
#LINKER		:= icx
#LINKER		:= gcc

.PHONY:	build debug build_subsystem debug_subsystem clean clean_build

build: 	build_subsystem
debug: 	debug_subsystem

clean_build: clean build_subsystem

build_subsystem:
	cd ./main		&& 		$(MAKE) build
	cd ./solver 	&& 		$(MAKE) build
	cd ./input 		&& 		$(MAKE) build

debug_subsystem:
	cd ./main		&& 		$(MAKE) debug
	cd ./solver 	&& 		$(MAKE) debug
	cd ./input 		&& 		$(MAKE) debug

clean:
	cd ./main		&& 		$(MAKE) clean
	cd ./solver 	&& 		$(MAKE) clean
	cd ./input 		&& 		$(MAKE) clean

install: 
	$(LINKER) $(LINKFLAGS) $(EXE_DIR)solver_ode ./*/*/*.o $(LINK_LIB)
	cp $(EXE_DIR)solver_ode /usr/local/bin/solver_ode_V3
	rm $(EXE_DIR)solver_ode


install_MKL:
	$(LINKER) $(LINKFLAGS) $(EXE_DIR)solver_ode ./*/*/*.o -lgsl -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lm -fopenmp
	cp $(EXE_DIR)solver_ode $(HOME)/bin/solver_ode_V3
	rm $(EXE_DIR)solver_ode

install_debug: debug
	$(LINKER) $(LINKFLAGS) $(EXE_DIR)solver_ode ./*/*/*.o $(LINK_LIB)
	cp $(EXE_DIR)solver_ode ./debug/solver_ode

clean_debug: clean
	rm $(EXE_DIR)solver_ode

uninstall:
	rm $(HOME)/bin/solver_ode_V3
