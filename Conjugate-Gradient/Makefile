CC = g++
CFLAGS = -std=c++11 -O2 -lm -c -Wall

all: conjGradient
	
conjGradient: conjGradient.o matrix.o main.o
	$(CC) -o conjGradient conjGradient.o matrix.o main.o

matrix.o: matrix.cpp
	$(CC) $(CFLAGS) matrix.cpp

conjGradient.o: conjGradient.cpp
	$(CC) $(CFLAGS) conjGradient.cpp

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

clean:
	rm -rf *o conjGradient
