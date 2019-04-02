# Directories
DIR=$(PWD)/
EXECDIR=$(DIR)exec/
OBJDIR=$(DIR)obj/
SRCDIR=$(DIR)src/

# GNU C compiler 
CC=mpicc

# Compiler flags
CFLAGS= -O2 -fopenmp -Wall
LIBFLAGS = -lm -lgsl -lgslcblas

# Command definition
RM=rm -f

# sources for main
sources_main = $(SRCDIR)main.c

objects_main = BGK.o momentRoutines.o transportroutines.o poissonNonlinPeriodic.o gauss_legendre.o input.o io.o zBar.o initialize_sol.o mesh.o implicit.o TNB.o

pref_main_objects = $(addprefix $(OBJDIR), $(objects_main))

sources_postproc = $(SRCDIR)TNB_postprocess.c $(SRCDIR)TNB.c $(SRCDIR)io.c

# linking step
MultiBGK: $(pref_main_objects) $(sources_main)
	@echo "Building Multispecies BGK code"
	$(CC) $(CFLAGS) -o $(EXECDIR)MultiBGK_ $(sources_main) $(pref_main_objects) $(LIBFLAGS)

postProc: $(sources_postproc)
	@echo "Building TNB postprocessor"
	$(CC) $(CFLAGS) -o $(EXECDIR)postProc_ $(sources_postproc) $(LIBFLAGS)


all: MultiBGK postProc

$(OBJDIR)BGK.o : $(SRCDIR)BGK.c
	@echo "Compiling  $< ... " ; \
	if [ -f  $@ ] ; then \
		rm $@ ;\
	fi ; \
	$(CC)  -c $(CFLAGS)  $< -o $@ 2>&1 ;

$(OBJDIR)momentRoutines.o : $(SRCDIR)momentRoutines.c
	@echo "Compiling  $< ... " ; \
	if [ -f  $@ ] ; then \
		rm $@ ;\
	fi ; \
	$(CC)  -c $(CFLAGS)  $< -o $@ 2>&1 ;


$(OBJDIR)zBar.o : $(SRCDIR)zBar.c
	@echo "Compiling  $< ... " ; \
	if [ -f  $@ ] ; then \
		rm $@ ;\
	fi ; \
	$(CC)  -c $(CFLAGS)  $< -o $@ 2>&1 ;



$(OBJDIR)transportroutines.o : $(SRCDIR)transportroutines.c
	@echo "Compiling  $< ... " ; \
	if [ -f  $@ ] ; then \
		rm $@ ;\
	fi ; \
	$(CC)  -c $(CFLAGS)  $< -o $@ 2>&1 ;

$(OBJDIR)gauss_legendre.o : $(SRCDIR)gauss_legendre.c
	@echo "Compiling  $< ... " ; \
	if [ -f  $@ ] ; then \
		rm $@ ;\
	fi ; \
	$(CC)  -c $(CFLAGS)  $< -o $@ 2>&1 ;

$(OBJDIR)poissonNonlinPeriodic.o : $(SRCDIR)poissonNonlinPeriodic.c
	@echo "Compiling  $< ... " ; \
	if [ -f  $@ ] ; then \
		rm $@ ;\
	fi ; \
	$(CC)  -c $(CFLAGS)  $< -o $@ 2>&1 ;

$(OBJDIR)input.o : $(SRCDIR)input.c
	@echo "Compiling  $< ... " ; \
	if [ -f  $@ ] ; then \
		rm $@ ;\
	fi ; \
	$(CC)  -c $(CFLAGS)  $< -o $@ 2>&1 ;

$(OBJDIR)io.o : $(SRCDIR)io.c
	@echo "Compiling  $< ... " ; \
	if [ -f  $@ ] ; then \
		rm $@ ;\
	fi ; \
	$(CC)  -c $(CFLAGS)  $< -o $@ 2>&1 ;

$(OBJDIR)initialize_sol.o : $(SRCDIR)initialize_sol.c
	@echo "Compiling  $< ... " ; \
	if [ -f  $@ ] ; then \
		rm $@ ;\
	fi ; \
	$(CC)  -c $(CFLAGS)  $< -o $@ 2>&1 ;

$(OBJDIR)mesh.o : $(SRCDIR)mesh.c
	@echo "Compiling  $< ... " ; \
	if [ -f  $@ ] ; then \
		rm $@ ;\
	fi ; \
	$(CC)  -c $(CFLAGS)  $< -o $@ 2>&1 ;

$(OBJDIR)implicit.o : $(SRCDIR)implicit.c
	@echo "Compiling  $< ... " ; \
	if [ -f  $@ ] ; then \
		rm $@ ;\
	fi ; \
	$(CC)  -c $(CFLAGS)  $< -o $@ 2>&1 ;

$(OBJDIR)TNB.o : $(SRCDIR)TNB.c
	@echo "Compiling  $< ... " ; \
	if [ -f  $@ ] ; then \
		rm $@ ;\
	fi ; \
	$(CC)  -c $(CFLAGS)  $< -o $@ 2>&1 ;

clean:
	$(RM) $(OBJDIR)*.o 
	$(RM) $(EXECDIR)*_

