CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall -DVERBOSE
BIN=../bin

all: JOSEE

JOSEE: JOSEE.c
	$(CC) $(CFLAGS) commonFunctions.c -lm comparisonFunctions.c evolutionaryEventsFunctions.c JOSEE.c -lpthread -o $(BIN)/JOSEE

clean:
	rm -rf $(BIN)/JOSEE
