#!/bin/bash
PNAME="./bin/Assignment4"

if [ "$#" -ne 3 ]; then
	echo "[ERROR] Illegal number of parameters"
	exit 1
fi

if [ ! -f $PNAME ]; then
	./rebuild_project.sh
fi

valgrind --leak-check=full ./$PNAME $1 $2 $3
