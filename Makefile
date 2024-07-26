COMPILATION = "C"

# Compiler
CC = gcc

# Folder structure
BIN		:= ./bin/
DEST	:= ./obj/
C		:= ./src/c/
OMP		:= ./src/omp/
CUDA	:= ./src/cuda/

# Output executable
TARGET = $(BIN)wildfires

# Commands
RM     = rm
MKDIR  = mkdir -p

# Flags
CFLAGS		= -Wall -std=c99 -lm -lfftw3 -O3
OMPFLAGS	= -fopenmp -lgomp
FLAGS		= $(CFLAGS) #$(OMPFLAGS)

# Source files
SRC_C = $(wildcard $(C)*.c)
SRC_OMP = $(wildcard $(OMP)*.c)
SRCS = $(SRC_C)
ifeq ($(COMPILATION), "OMP")
# Exclude files in SRC_OMP from SRC
EXCLUDE_SRC = $(C)pde.c $(C)parameters.c
SRCS = $(filter-out $(EXCLUDE_SRC),$(SRC_C))
SRCS += $(SRC_OMP)
FLAGS += $(OMPFLAGS)
endif
$(info SRCS = $(SRCS))

# Object files
OBJS = $(patsubst %.c,$(DEST)%.o,$(notdir $(SRCS)))
$(info OBJS = $(OBJS))

# Default target
all: directories $(TARGET)

# Link the executable
$(TARGET): $(OBJS)
	$(CC) -o $@ $(OBJS) $(FLAGS)

# Compile source files into object files
$(DEST)%.o: $(C)%.c
	$(CC) $(FLAGS) -c $< -o $@

# # Object files
# OBJS = $(SRCS:.c=.o)

# # Default target
# all: directories $(TARGET)

# # Link the executable
# $(TARGET): $(OBJS)
# 	$(CC) -o $@ $(OBJS) $(FLAGS)

# # Compile source files into object files
# %.o: %.c
# 	$(CC) $(FLAGS) -c $< -o $@

directories:
	$(MKDIR) $(BIN)
	$(MKDIR) $(DEST)

# # Clean up
clean:
	$(RM) -rf $(DEST)*.o

distclean: clean
	$(RM) -rf $(BIN)*

.PHONY: all clean
