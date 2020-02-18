########################################################################
####################### Makefile Template ##############################
########################################################################

include ../../make.inc
# compiler flags
FFLAGS = $(F90FLAGS) $(OMP)
LDFLAGS =

# Makefile settings - Can be customized.
APPNAME = sym.a
EXT = .F90
SRCDIR = ./
OBJDIR = obj

MODULES= modelk.F90
SRC = $(MODULES) findsymcrys.F90 findsym.F90 findsymlat.F90  symmetry.F90  writesym.F90


############## Do not change anything from here downwards! #############
OBJ = $(SRC:%.F90=$(OBJDIR)/%.o)
MOD = $(MODULES:%.F90=%.mod)
RM = rm
DELOBJ = $(OBJ)

########################################################################
####################### Targets beginning here #########################
########################################################################

all: $(APPNAME)

# Builds the app
$(APPNAME): $(OBJ)
	ar -rc $@ $^ $(LDFLAGS)

# Building rule for .o files and its .c/.cpp in combination with all .h
$(OBJDIR)/%.o: $(SRCDIR)/%$(EXT)
	$(F90) $(FFLAGS) -o $@ -c $<

################### Cleaning rules for Unix-based OS ###################
# Cleans complete project
.PHONY: clean
clean:
	$(RM) $(DELOBJ) $(DEP) $(APPNAME) $(MOD)

