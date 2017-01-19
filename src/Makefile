CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall -DVERBOSE
BIN=../bin

all: JOSEE cutter

JOSEE: JOSEE.cpp
	$(CXX) $(CFLAGS) commonFunctions.cpp -lm comparisonFunctions.cpp evolutionaryEventsFunctions.cpp class_structs.cpp JOSEE.cpp -lpthread -o $(BIN)/JOSEE

cutter: cutter.cpp
	$(CXX) $(CFLAGS) commonFunctions.cpp -lm comparisonFunctions.cpp cutter.cpp -lpthread -o $(BIN)/cutter

clean:
	rm -rf $(BIN)/JOSEE $(BIN)/cutter
