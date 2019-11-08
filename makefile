# Directories
DIR=$(PWD)/
EXECDIR=$(DIR)exec/
OBJDIR=$(DIR)obj/
SRCDIR=$(DIR)src/

# GNU C compiler 
CC=mpicc

# Compiler flags
CFLAGS= -O0  -Wall -g

LIBFLAGS = -lm -lgsl -lgslcblas -fopenmp

# Command definition
RM=rm -f

# sources for main
sources_main = $(SRCDIR)main.c


objects_main = BGK.o momentRoutines.o transportroutines.o poissonNonlinPeriodic.o gauss_legendre.o input.o io.o zBar.o initialize_sol.o mesh.o implicit.o poissonNonlinNonPeriodic.o parallel_poison.o


pref_main_objects = $(addprefix $(OBJDIR), $(objects_main))


# linking step
MultiBGK: $(pref_main_objects) $(sources_main)
	@echo "Building Multispecies BGK code"
	$(CC) $(CFLAGS) -o $(EXECDIR)MultiBGK_ $(sources_main) $(pref_main_objects) $(LIBFLAGS)



$(OBJDIR)%.o : $(SRCDIR)%.c
	@echo "Compiling  $< ... " ; \
	if [ -f  $@ ] ; then \
		rm $@ ;\
	fi ; \
	$(CC)  -c $(CFLAGS)  $< -o $@ 2>&1 ;

clean:
	$(RM) $(OBJDIR)*.o 
	$(RM) $(EXECDIR)*_

