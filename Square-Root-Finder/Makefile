CC = g++
CFLAGS = -c -Wall

all: main

main: main.o raiz.o
	$(CC) main.o raiz.o -o main

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

raiz.o: raiz.cpp
	$(CC) $(CFLAGS) raiz.cpp

clean:
	rm -rf *o main
