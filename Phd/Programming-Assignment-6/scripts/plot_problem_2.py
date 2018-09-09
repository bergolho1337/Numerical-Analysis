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

def showPoints (x,y,colorname,labelname):
    #plt.ylim(0,10)
    plt.scatter(x,y,label=labelname,marker='o',s=40,c=colorname)

def plotSolution (x,y,colorname,labelname):
    #plt.ylim(0,10)
    plt.plot(x,y,label=labelname,c=colorname)

def showOutput ():
    plt.grid()
    plt.xlabel(u"x",fontsize=15)
    plt.ylabel(u"y",fontsize=15)
    plt.title(u"Curva ajustada por Mínimos Quadrados (Problema 2)",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.show()

def writeOutput ():
    plt.grid()
    plt.xlabel(u"x",fontsize=15)
    plt.ylabel(u"y",fontsize=15)
    plt.title(u"Curva ajustada por Mínimos Quadrados (Problema 2)",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.savefig("output/problem_2.pdf")

def plotLimitPoint (x,y):
    plt.scatter(x,y,label="limit",marker='x',s=60,c="red")

def main():

    points_filename = sys.argv[1]
    points_filename2 = sys.argv[2]
    solution_filename = sys.argv[3]
    solution_filename2 = sys.argv[4]
    
    xpts1, ypts1 = readPoints(points_filename)
    xpts2, ypts2 = readPoints(points_filename2)
    
    x1, y1 = readPoints(solution_filename)
    x2, y2 = readPoints(solution_filename2)
    
    showPoints(xpts1,ypts1,"red","set_1")
    plotSolution(x1,y1,"red","set_1")

    showPoints(xpts2,ypts2,"blue","set_2")
    plotSolution(x2,y2,"blue","set_2")

    writeOutput()
    #showOutput()
    #plt.show()

if __name__ == "__main__":
    main()
