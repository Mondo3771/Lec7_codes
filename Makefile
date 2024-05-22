INC="./inc"
FLAGS=-I$(INC)
OMPFLAGS=-fopenmp
CFLAGS = -g -Wall
CC=gcc
MPICC=mpicc
CPP=g++
	
#SOURCES=$(wildcard *.c)
OBJECTS=$(patsubst %.c, %, $(wildcard *.c))


all: $(OBJECTS)
	@echo 'objects are "$(OBJECTS)"'

%: %.c
	$(MPICC) $(CFLAGS) $< -o $@ -lm
	
clean:
	rm -vf $(OBJECTS)
