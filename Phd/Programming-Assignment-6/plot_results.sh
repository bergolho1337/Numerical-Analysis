#!/bin/bash

if [ "$#" -ne 1 ]; then
	echo "***********************************************"
    echo "Usage:> ./plot_results.sh <problem_id>"
    echo "***********************************************"
    echo "<problem_id> = Problem identifier"
    echo "            1 - Problem 1"
    echo "            2 - Problem 2"
    echo "***********************************************"
    exit 1
fi

if [ ! -f $PNAME ]; then
	./rebuild_project.sh
fi

if [ $1 -eq 1 ]; then
    echo "[Problem 1] Generating plot ..."
	python ./scripts/plot_problem_1.py ./inputs/points-problem-1.txt ./solution.dat
    echo "[Problem 1] Done"
elif [ $1 -eq 2 ]; then
    echo "[Problem 2] Generating plot ..."
	 python ./scripts/plot_problem_2.py ./inputs/points1-problem-2.txt ./inputs/points2-problem-2.txt ./solution.dat ./solution2.dat
    echo "[Problem 2] Done"
fi


