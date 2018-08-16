# =====================================================================
# Script to convert the Heder files to the MatrixMarket format
# =====================================================================

import sys
import numpy as np

def readMatrix (filename):
    A = np.loadtxt(filename)
    return A

def readRHS (filename):
    b = np.loadtxt(filename)
    return b

def printMatrix (A):
    rows = A.shape[0]
    cols = A.shape[1]

    print("A:\n")
    for i in range(rows):
        for j in range(cols):
            print("%.10lf " % A[i][j]),
        print("\n")

def printRHS (b):
    size = b.shape[0]

    print("b:\n")
    for i in range(size):
        print("%.10lf\n" % b[i])

def writeMatrix (A):
    rows = A.shape[0]
    cols = A.shape[1]

    file = open("output/A.mtx","w")
    file.write("%%MatrixMarket matrix coordinate real general\n")
    file.write("%d %d %d\n" % (rows,cols,rows*cols))
    for i in range(rows):
        for j in range(cols):
            file.write("%d %d %g\n" % (i+1,j+1,A[i][j]))
    file.close()
    print("[+] Matrix saved sucessfully in file:> \"output/A.mtx\"")

def writeRHS (b):
    size = b.shape[0]

    file = open("output/b.mtx","w")
    file.write("%%MatrixMarket matrix array real general\n")
    file.write("%d %d\n" % (size,1))
    for i in range(size):
        file.write("%g\n" % (b[i]))
    print("[+] RHS saved sucessfully in file:> \"output/b.mtx\"")

if (len(sys.argv) != 3):
    print("----------------------------------------")
    print("Usage:> python converter.py <A> <b>")
    print("----------------------------------------")
    print("<A> = Coefficient matrix filename")
    print("<b> = Right-Hand-Side array filename")
    print("----------------------------------------")
else:
    matrix_filename = sys.argv[1]
    rhs_filename = sys.argv[2]

    A = readMatrix(matrix_filename)
    b = readRHS(rhs_filename)

    #printMatrix(A)
    #printRHS(b)

    writeMatrix(A)
    writeRHS(b)
