CC = gcc
CFLAGS = -c -Wall

all: main

main: main.o
	$(CC) main.o -o main -lm

main.o: main.c
	$(CC) $(CFLAGS) main.c -lm

clean:
	rm -rf *o main
