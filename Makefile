CC = g++
CFLAGS = -Wall -g

all: gillespie

gillespie: gillespie.cpp
	$(CC) $(CFLAGS) -o $(@) $(^)
