CC = g++
CFLAGS = -c -Wall

all: euler

euler: euler.o 
	$(CC) euler.o -o euler

euler.o: euler.cpp
	$(CC) $(CFLAGS) euler.cpp

clean:
	rm -rf *o euler
