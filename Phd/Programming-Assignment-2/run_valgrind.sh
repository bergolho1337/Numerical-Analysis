#!/bin/bash
PNAME="./bin/Assignment2"

if [ "$#" -ne 3 ]; then
	echo "[ERROR] Illegal number of parameters"
	exit 1
fi

if [ ! -f $PNAME ]; then
	./recompile_project.sh
fi

# Run valgrind to check for heap memory leaks ...
valgrind --leak-check=full --show-leak-kinds=all ./$PNAME $1 $2 $3
