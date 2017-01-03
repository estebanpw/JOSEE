CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall -DVERBOSE
BIN=../bin

all: JOSEE

JOSEE: JOSEE.cpp
	$(CXX) $(CFLAGS) commonFunctions.cpp -lm comparisonFunctions.cpp evolutionaryEventsFunctions.cpp class_structs.cpp JOSEE.cpp -lpthread -o $(BIN)/JOSEE

clean:
	rm -rf $(BIN)/JOSEE
