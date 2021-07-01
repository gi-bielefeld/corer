CC = g++
CFLAGS = -march=native -lbifrost -pthread -lz -lpthread -std=c++11 -O3

all:
	$(CC) Corer.cpp $(CFLAGS) -o Corer
