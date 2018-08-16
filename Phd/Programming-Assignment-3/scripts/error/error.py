# =====================================================================
# Program that computes the error using the infinite norm
# =====================================================================

import sys
import numpy as np

# Read the solution file
def readSolution (filename):
    s = np.loadtxt(filename)
    return s

# Answer array from Heder's handout
def compAnswer (s):
    size = len(s)
    ans = np.zeros(size)
    for i in range(size):
        ans[i] = 101.0
    return ans

# Calculates the norm-inf of between the answer and the solution array
def calcErrorInf (s,a):
    size = len(s)
    norm = abs(s[0]-a[0])
    for i in range(size-1):
        value = abs(s[i+1]-a[i+1])
        if (value > norm):
            norm = value
    return norm



if (len(sys.argv) != 2):
    print("----------------------------------------")
    print("Usage:> python error.py <solution>")
    print("----------------------------------------")
    print("<solution> = Solution filename")
    print("----------------------------------------")
else:
    solution_filename = sys.argv[1]

    sol = readSolution(solution_filename)
    ans = compAnswer(sol)
    error = calcErrorInf(sol,ans)

    print("****************************************")
    print("Error = %g" % error)
    print("****************************************")
    #for i in range(len(sol)):
    #    print("%g\n" % (abs(sol[i]-ans[i])))
