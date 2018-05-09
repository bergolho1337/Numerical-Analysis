CC = g++
CFLAGS = -c -Wall

all: main

main: main.o sistema.o matriz.o
	$(CC) main.o sistema.o matriz.o -o main

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

sistema.o: sistema.cpp
	$(CC) $(CFLAGS) sistema.cpp

matriz.o: matriz.cpp
	$(CC) $(CFLAGS) matriz.cpp

clean:
	rm -rf *o main
