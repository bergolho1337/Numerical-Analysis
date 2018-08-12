# This program calculates the forward and backward error of an aproximation of the sine function considering only the first term of its 
# corresponding series, which is given by the formula:
#
# 	sin(x) = x - (x^3 / 3!) + (x^5 / 5!) - (x^7 / 7!) + ...
#
# y = sin(x) || y^ = x
#
# This program also solves the problem 1.6 from Michael Heath Scientific Computing book. 

import sys
import math

def calc_errors (x,y,nterm):
	print("[!] Using %d term of the series ..." % nterm)
	if (nterm == 1):
		y_aprox = x
	elif (nterm == 2):
		y_aprox = x - (math.pow(x,3) / 6.0)
	else:
		print("Error !")		
		
	x_aprox = math.asin(y_aprox)

	# Forward error
	fe = y_aprox - y	

	# Backward error
	be = x_aprox - x

	print("y_aprox = %.10lf" % y_aprox)
	print("x_aprox = %.10lf" % x_aprox)

	return fe, be

if (len(sys.argv) != 3):
	print("Usage:> python calculator.py <x> <nterm>")
else:
	x = float(sys.argv[1])
	nterm = int(sys.argv[2])

	y = math.sin(x)

	fe, be = calc_errors(x,y,nterm)		

	print("x = %.10lf" % x)
	print("y = %.10lf" % y)

	print("Forward error = %g" % fe)
	print("Backward error = %g" % be)


