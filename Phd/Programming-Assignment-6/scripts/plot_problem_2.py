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

def readSolution (filename):
    file = open(filename,"r")

    line = file.readline()
    n = int(line)

    x = np.zeros(n)
    y = np.zeros(n)
    z = np.zeros(n)

    for i in range(n):
        line = file.readline()
        tokens = line.split()

        x[i] = float(tokens[0])
        y[i] = float(tokens[1])
        z[i] = float(tokens[2])

    return x, y, z

def showPoints (x,y,colorname,labelname):
    #plt.ylim(0,10)
    plt.scatter(x,y,label=labelname,marker='o',s=40,c=colorname)

def plotSolution (x,y,z,colorname1,colorname2,labelname1,labelname2):
    #plt.ylim(0,10)
    plt.plot(x,y,label=labelname1,c=colorname1)
    plt.plot(x,z,label=labelname2,c=colorname2)

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
    
    xpts1, ypts1 = readPoints(points_filename)
    xpts2, ypts2 = readPoints(points_filename2)
    
    x1, y1, z1 = readSolution(solution_filename)
    
    showPoints(xpts1,ypts1,"red","set_1")
    showPoints(xpts2,ypts2,"blue","set_2")
    plotSolution(x1,y1,z1,"red","blue","set_1","set_2")

    writeOutput()
    #showOutput()
    #plt.show()

if __name__ == "__main__":
    main()
