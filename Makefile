COMPILATION = "CUDA"

# Compiler
CC = g++
NVCC = nvcc

# Folder structure
SRC     := ./src/
BIN		:= ./bin/
DEST	:= ./obj/
CSRC	:= $(SRC)c/

# Output executable
EXE	= wildfires
TARGET	= $(BIN)$(EXE)

# Commands
RM     = rm
MKDIR  = mkdir -p

# Flags
CFLAGS		= -Wall -lm -lfftw3 -O3 #-std=c99
OMPFLAGS	= -fopenmp -lgomp
NVARCH		= 75
NVFLAGS		= -arch=sm_$(NVARCH) -Xptxas -dlcm=ca
FLAGS		= $(CFLAGS)

# Source files
MAIN   = cu/main.cu 
CODC = output.c parameters.c logs.c
# CODC = $(wildcard $(CSRC)*.c) 
CODOMP = #$(wildcard $(OMP)*.c)
CODCU = functions.cu pde.cu solver.cu turbulence.cu utils.cu

# Folders
FC     = c/
FCU    = cu/

# Object files
# OBJS = $(patsubst %.c,$(DEST)%.o,$(notdir $(SRCS)))
OBJC   = $(patsubst %.c,$(DEST)$(FC)%.o,$(notdir $(CODC)))
OBJCU  = $(patsubst %.cu,$(DEST)$(FCU)%.o,$(notdir $(CODCU)))
# $(info $(OBJC))
# $(info $(OBJCU))

SRCMAIN = $(patsubst %,$(SRC)%,$(MAIN))
OBJMAIN = $(patsubst $(SRC)%.cu,$(DEST)%.o,$(SRCMAIN))
# $(info $(SRCMAIN))
# $(info $(OBJMAIN))

# Default target
all: directories $(BIN)$(EXE)

.PHONY: clean

$(BIN)$(EXE): $(OBJC) $(OBJCU) $(OBJMAIN)
	$(NVCC) $(NVFLAGS) $^ -o $@

# $(BIN)$(EXE): $(OBJC) $(OBJCU) $(OBJMAIN)
# 	$(CC) $^ -o $@ $(CFLAGS)

$(OBJMAIN): $(SRCMAIN)
	$(NVCC) $(NVFLAGS) -c $? -o $@

$(OBJC): $(DEST)%.o : $(SRC)%.c
	$(CC) $(CFLAGS) -c $? -o $@

$(OBJCU): $(DEST)%.o : $(SRC)%.cu
	$(NVCC) $(NVFLAGS) -c $? -o $@

directories:
	$(MKDIR) $(BIN)
	$(MKDIR) $(DEST)$(FC) 
	$(MKDIR) $(DEST)$(FCU)  

# # Clean up
clean:
	$(RM) -rf $(DEST)$(FC)*.o
	$(RM) -rf $(DEST)$(FCU)*.o

distclean: clean
	$(RM) -rf $(BIN)*
