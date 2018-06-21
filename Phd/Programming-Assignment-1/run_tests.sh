#!/bin/bash
PROGRAM_NAME="assignment1"
XVALUES=( 1 2 10 20 -1 -2 -10 -20 )

if [ ! -f $PROGRAM_NAME ]; then
    make
else
    make remade
fi

for X in ${XVALUES[@]}; do
    echo "------------------------------------------------------"
    echo "Running test for x = $X"
    ./$PROGRAM_NAME $X
    echo "------------------------------------------------------"
done

