########################################################################
####################### Makefile Template ##############################
########################################################################

include ../make.inc

# additional libraries
SYMLIB=./SYM/sym.a
ARPACKLIB=./ARPACK/libarpack.a 

# compiler flags
FFLAGS = $(F90FLAGS) $(OMP) $(DPARDI) $(DMPI)

# linkink library flags
LDFLAGS = $(SYMLIB) $(ARPACKLIB) $(BLAS) $(LAPACK) 

# modules
MODULES = modcom.F90 geometry_library.F90 parameters.F90 symmetryclass.F90 gridclass.F90  \
          slater_koster.F90 tbclass.F90 tasksclass.F90
# sourses
SRC     = $(MODULES) arpack_interface.F90  main.F90

# Makefile settings
APPNAME = tbx
OBJDIR  = obj
EXT     = .F90

############## Do not change anything from here downwards! #############
OBJ = $(SRC:%.F90=$(OBJDIR)/%.o)
MOD = $(MODULES:%.F90=%.mod)
RM = rm
DELOBJ = $(OBJ)

########################################################################
####################### Targets beginning here #########################
########################################################################


all: arpack sym $(APPNAME) 

arpack: 
	cd ARPACK; make lib; cd ..;
sym: 
	cd SYM; make; cd ..;

# Builds the app
$(APPNAME): $(OBJ)
	$(F90) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Building rule for .o files and its .c/.cpp in combination with all .h
$(OBJDIR)/%.o: %$(EXT)
	$(F90) $(FFLAGS) -o $@ -c $<

clean:
	cd SYM; make clean; cd ..;
	$(RM) $(DELOBJ) $(DEP) $(APPNAME) $(MOD)
cleanall:
	cd ARPACK;  make clean; cd ..;
	cd SYM; make clean; cd ..;
	$(RM) $(DELOBJ) $(DEP) $(APPNAME) $(MOD)

# create docs
docs:
	doxygen Doxyfile
