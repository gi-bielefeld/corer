CC = g++
CFLAGS = -DMAX_KMER_SIZE=64 -march=native -lbifrost -pthread -lz -lpthread -std=c++11 -Wall -O3

all:
	$(CC) Corer.cpp $(CFLAGS) -o Corer
