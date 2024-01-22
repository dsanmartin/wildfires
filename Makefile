#
# Folder structure
#

BIN     := ./bin/
SRC     := ./src/
DEST    := ./obj/

EXE    = wildfires

#
# Executables
#

CC     = gcc
RM     = rm
MKDIR  = mkdir -p

#
# C/C++ flags
#

CFLAGS    = -Wall -g -O3 -std=c99 -lm -lfftw3#-Wextra

#
# Files to compile: 
#

MAIN   = main.c
CODC   = functions.c logs.c output.c parameters.c pde.c poisson.c turbulence.c utils.c

#
# Formating the folder structure for compiling/linking/cleaning.
#

FC     = 

#
# Preparing variables for automated prerequisites
#

OBJC   = $(patsubst %.c,$(DEST)$(FC)%.o,$(CODC))

SRCMAIN = $(patsubst %,$(SRC)%,$(MAIN))
OBJMAIN = $(patsubst $(SRC)%.c,$(DEST)%.o,$(SRCMAIN))

# .PHONY: directories

#
# The MAGIC
#
 all: directories $(BIN)$(EXE)

.PHONY: clean

$(BIN)$(EXE): $(OBJC) $(OBJMAIN)
	$(CC) $^ -o $@ $(CFLAGS)

$(OBJMAIN): $(SRCMAIN)
	$(CC) $(CFLAGS) -c $? -o $@

$(OBJC): $(DEST)%.o : $(SRC)%.c
	$(CC) $(CFLAGS) -c $? -o $@


directories:
	$(MKDIR) $(BIN)
	$(MKDIR) $(DEST)$(FC) 

#
# Makefile for cleaning
# 

clean:
	$(RM) -rf $(DEST)*.o

fresh:
	$(RM) -rf test/output/*

distclean: clean
	$(RM) -rf $(BIN)*


