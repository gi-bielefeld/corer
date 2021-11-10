CC = g++
CFLAGS = -march=native -std=c++11 -O3
LINKS = -lbifrost -pthread -lz -lpthread -lgtest -lgtest_main
SRCFILES = gtest_main.cpp IO_unittest.cpp Core_unittest.cpp Bridging_unittest.cpp Traversal_unittest.cpp

all:
	$(CC) $(CFLAGS) -o MyTest $(SRCFILES) $(LINKS) -I../src
