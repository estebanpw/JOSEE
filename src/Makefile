CC=gcc
CXX=g++ -std=c++11
CFLAGS=-g -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall -Wno-literal-suffix -DVERBOSE
BIN=../bin

all: JOSEE cutter

JOSEE: JOSEE.cpp
	$(CXX) $(CFLAGS) alignment_functions.cpp -lm commonFunctions.cpp -lm comparisonFunctions.cpp evolutionaryEventsFunctions.cpp class_structs.cpp JOSEE.cpp -lpthread -o $(BIN)/JOSEE

cutter: cutter.cpp
	$(CXX) $(CFLAGS) alignment_functions.cpp -lm commonFunctions.cpp -lm comparisonFunctions.cpp class_structs.cpp cutter.cpp -lpthread -o $(BIN)/cutter

evengen: event_generator.cpp
	$(CXX) $(CFLAGS) alignment_functions.cpp -lm commonFunctions.cpp -lm comparisonFunctions.cpp class_structs.cpp evolution.cpp event_generator.cpp -lpthread -o $(BIN)/evengen


clean:
	rm -rf $(BIN)/JOSEE $(BIN)/cutter $(BIN)/evengen
