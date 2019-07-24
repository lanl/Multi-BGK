# Directories
DIR=$(PWD)/
EXECDIR=$(DIR)exec/
OBJDIR=$(DIR)obj/
SRCDIR=$(DIR)src/
ALDR_GLUE_DIR=/Users/haack/ALDR/glueCode_proofOfConcept/
SQL_DIR=/Users/haack/spack/opt/spack/darwin-highsierra-x86_64/clang-4/sqlite-3.26.0-pbpvhmzw5xp5abrxclqp4tjmdbeht2dz
# GNU C compiler 
CC=mpicc

# Compiler flags
CFLAGS= -O0 -fopenmp -Wall
LIBFLAGS = -lm -lgsl -lgslcblas

# Command definition
RM=rm -f

# sources for main
sources_main = $(SRCDIR)main.c

sources_aldr = $(sources_main) $(ALDR_GLUE_DIR)alInterface.cpp
#sources_aldr = $(sources_main)

objects_main = BGK.o momentRoutines.o transportroutines.o poissonNonlinPeriodic.o gauss_legendre.o input.o io.o zBar.o initialize_sol.o mesh.o implicit.o

pref_main_objects = $(addprefix $(OBJDIR), $(objects_main))

# linking step
MultiBGK: $(pref_main_objects) $(sources_main)
	@echo "Building Multispecies BGK code"
	$(CC) $(CFLAGS) -o $(EXECDIR)MultiBGK_ $(sources_main) $(pref_main_objects) $(LIBFLAGS)

ALDR: CFLAGS += -DALDR_ON -I$(ALDR_GLUE_DIR) -DSQLITE_DIR=$(SQL_DIR)

ALDR: $(pref_main_objects) $(sources_aldr) 
	@echo "Building Multispecies BGK code"
	$(CC) $(CFLAGS) -o $(EXECDIR)MultiBGK_AL_ $(sources_aldr) $(pref_main_objects) $(LIBFLAGS)

$(OBJDIR)%.o : $(SRCDIR)%.c
	@echo "Compiling  $< ... " ; \
	if [ -f  $@ ] ; then \
		rm $@ ;\
	fi ; \
	$(CC)  -c $(CFLAGS)  $< -o $@ 2>&1 ;



clean:
	$(RM) $(OBJDIR)*.o 
	$(RM) $(EXECDIR)*_