#!/bin/bash
PNAME="./bin/Assignment6"

if [ "$#" -ne 5 ]; then
	echo "[ERROR] Illegal number of parameters"
	exit 1
fi

if [ ! -f $PNAME ]; then
	./rebuild_project.sh
fi

valgrind --leak-check=full --show-leak-kinds=all ./$PNAME $1 $2 $3 $4 $5
