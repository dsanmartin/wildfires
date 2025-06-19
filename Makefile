COMPILATION = "CUDA"

# Compiler
CC = gcc
NVCC = nvcc

# Folder structure
SRC     := ./src/
BIN		:= ./bin/
DEST	:= ./obj/

# Output executable
EXE	= wildfires

# Commands
RM     = rm
MKDIR  = mkdir -p

# Flags
OMPFLAGS	= -fopenmp -lgomp
NVARCH		= 89
NVFLAGS		= -gencode arch=compute_$(NVARCH),code=sm_$(NVARCH) -lcufft -O3# -lcufft -Xptxas -dlcm=ca -I/usr/local/cuda/inc -L/usr/local/cuda/lib
CFLAGS		= -Wall -lm -O3 -lfftw3

# Source files
ifeq ($(COMPILATION), "CUDA")
EXE = wildfires_cu
MAIN   = cu/main.cu 
CODC = #poisson2.c
# CODC = $(wildcard $(CSRC)*.c) 
CODOMP = #$(wildcard $(OMP)*.c)
CODCU = functions.cu ibm.cu logs.cu input.cu output.cu parameters.cu pde.cu pressure.cu solver.cu turbulence.cu utils.cu
else ifeq ($(COMPILATION), "C")
# CFLAGS	+= -lfftw3
EXE = wildfires_c
MAIN   = c/main.c
CODC = functions.c ibm.c logs.c output.c parameters.c pde.c poisson.c solver.c turbulence.c utils.c
endif

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

ifeq ($(COMPILATION), "CUDA")

$(BIN)$(EXE): $(OBJC) $(OBJCU) $(OBJMAIN)
	$(NVCC) $(NVFLAGS) $^ -o $@

$(OBJMAIN): $(SRCMAIN)
	$(NVCC) $(NVFLAGS) -c $? -o $@

$(OBJC): $(DEST)%.o : $(SRC)%.c
	$(CC) $(CFLAGS) -c $? -o $@

$(OBJCU): $(DEST)%.o : $(SRC)%.cu
	$(NVCC) $(NVFLAGS) -c $? -o $@

else ifeq ($(COMPILATION), "C")

$(BIN)$(EXE): $(OBJC) $(OBJMAIN)
	$(CC) $^ -o $@ $(CFLAGS)

# $(OBJMAIN): $(SRCMAIN)
# 	$(CC) $(CFLAGS) -c $? -o $@

$(OBJC): $(DEST)%.o : $(SRC)%.c
	$(CC) $(CFLAGS) -c $? -o $@

endif

directories:
	$(MKDIR) $(BIN)
	$(MKDIR) $(DEST)$(FC) 
	$(MKDIR) $(DEST)$(FCU)
	$(MKDIR) data/output/  

# # Clean up
clean:
	$(RM) -rf $(DEST)$(FC)*.o
	$(RM) -rf $(DEST)$(FCU)*.o

distclean: clean
	$(RM) -rf $(BIN)*
