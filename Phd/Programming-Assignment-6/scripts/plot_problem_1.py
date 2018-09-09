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

def showPoints (x,y):
    #plt.ylim(0,10)
    plt.scatter(x,y,label="dataset",marker='o',s=40)

def plotSolution (x,y):
    #plt.ylim(0,10)
    plt.plot(x,y,label="aprox")

def showOutput ():
    plt.grid()
    plt.xlabel(u"x",fontsize=15)
    plt.ylabel(u"y",fontsize=15)
    plt.title(u"Curva ajustada por Mínimos Quadrados (Problema 1)",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.show()

def writeOutput ():
    plt.grid()
    plt.xlabel(u"x",fontsize=15)
    plt.ylabel(u"y",fontsize=15)
    plt.title(u"Curva ajustada por Mínimos Quadrados (Problema 1)",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.savefig("problem_1.pdf")

def plotLimitPoint (x,y):
    plt.scatter(x,y,label="limit",marker='x',s=60,c="red")

def main():

    points_filename = sys.argv[1]
    solution_filename = sys.argv[2]
    
    x1, y1 = readPoints(points_filename)
    x2, y2 = readPoints(solution_filename)
    x3, y3 = [11.6160000000, 1998.8567991085]
    
    showPoints(x1,y1)
    plotSolution(x2,y2)
    plotLimitPoint(x3,y3)

    writeOutput()
    #showOutput()
    #plt.show()

if __name__ == "__main__":
    main()
