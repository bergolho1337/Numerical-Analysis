# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt

def readPoints (filename):
    file = open(filename,"r")

    line = file.readline()
    n = int(line)

    x = np.zeros(n)
    y = np.zeros(n)

    for i in range(n):
        line = file.readline()
        tokens = line.split()

        x[i] = float(tokens[0])
        y[i] = float(tokens[1])

    return x, y

def buildLine (x,y,i,j):
    a = (y[j] - y[i])/(x[j] - x[i])
    b = -a * x[i] + y[i]
    return a, b

def showPoints (x,y):
    #plt.ylim(0,10)
    plt.scatter(x,y,label="dataset",marker='o',s=40,color="black")

def showOutput ():
    plt.grid()
    plt.xlabel(u"x",fontsize=15)
    plt.ylabel(u"y",fontsize=15)
    plt.title(u"Curva ajustada por Mínimos Quadrados (Problema 3)",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.show()

def writeOutput ():
    plt.grid()
    plt.xlabel(u"x",fontsize=15)
    plt.ylabel(u"y",fontsize=15)
    plt.title(u"Curva ajustada (Problema 3)",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.savefig("../output/problem_3.pdf")

def showLine (x,y,i,j,labelname,colorname):
    a, b = buildLine(x,y,i,j)

    print("(%d,%d)" % (i,j))
    print("ax + b => %.10lf.x + %.10lf" % (a,b))
    X = np.linspace(0,8,100)
    Y = a*X + b

    plt.plot(X,Y,label=labelname,c=colorname)

def showSolutionLine (x,y,id1,id2,id3,id4,labelname,colorname):
    a1, b1 = buildLine(x,y,id1,id2)
    a2, b2 = buildLine(x,y,id3,id4)
    b3 = (b2 + b1) / 2.0

    error = evaluateSolution(x,y,a1,b3)

    print("Solution")
    print("ax + b => %.10lf.x + %.10lf" % (a1,b3))
    print("Error = %.10lf" % error)
    X = np.linspace(0,8,100)
    Y = a1*X + b3

    plt.plot(X,Y,label=labelname,c=colorname)

def printAllLines (x,y):
    n = len(x)
    for i in range(n):
        for j in range(i+1,n):
            if (i != j):
                a, b = buildLine(x,y,i,j)

                print("(%d,%d)" % (i,j))
                print("ax + b => %.10lf.x + %.10lf" % (a,b))
                print

def evaluateSolution (x,y,a,b):
    n = len(x)
    max_error = 0.0

    for i in range(n):
        delta_x = abs( (y[i] - b)/a - x[i] )
        delta_y = abs( a*x[i] + b - y[i] )
        error = max(delta_x,delta_y)
        if (error > max_error):
            max_error = error

    return max_error

def main():

    points_filename = sys.argv[1]
    
    x1, y1 = readPoints(points_filename)
    
    showPoints(x1,y1)
    
    # Debug purposes ...
    #printAllLines(x1,y1)
    
    # Solution algorithm
    showLine(x1,y1,1,3,"1,3","red")
    showLine(x1,y1,2,4,"2,4","blue")
    #showSolutionLine(x1,y1,0,2,1,4,"aprox","green")
    showSolutionLine(x1,y1,1,3,2,4,"aprox","green")

    writeOutput()
    #showOutput()

if __name__ == "__main__":
    main()
