BIN  = ./bin
SRC  = ./src
OBJ  = ./obj
EXE  = wildfires

CC = gcc
CFLAGS = -Wall -Wextra

all: project

project: $(OBJ)/main.o $(OBJ)/parameters.o $(OBJ)/logs.o
	$(CC) $(CFLAGS) -o $(BIN)/$(EXE) $(OBJ)/main.o $(OBJ)/parameters.o $(OBJ)/logs.o

$(OBJ)/main.o: $(SRC)/main.c
	$(CC) $(CFLAGS) -c $(SRC)/main.c -o $(OBJ)/main.o

$(OBJ)/parameters.o: $(SRC)/parameters.c
	$(CC) $(CFLAGS) -c $(SRC)/parameters.c -o $(OBJ)/parameters.o

$(OBJ)/logs.o: $(SRC)/logs.c
	$(CC) $(CFLAGS) -c $(SRC)/logs.c -o $(OBJ)/logs.o

clean:
	rm -f project $(OBJ)/main.o $(OBJ)/parameters.o $(OBJ)/logs.o
