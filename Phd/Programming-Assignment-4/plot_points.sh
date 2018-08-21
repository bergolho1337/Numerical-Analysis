#!/bin/bash
INPUT_FILE="../inputs/snoopy_points.txt"

cd scripts 
python plot.py $INPUT_FILE interpolate_points.txt