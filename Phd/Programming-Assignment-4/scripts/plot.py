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

def readInterpolatePoints (filename):
    z = np.loadtxt(filename)
    return z

def showPoints (x,y):
    #plt.ylim(0,10)
    plt.scatter(x,y)

def showInterpolatePoints (z, colorname, labelname):
    #plt.ylim(0,20)
    plt.plot(z[:,0],z[:,1],c=colorname,label=labelname)
    #plt.grid() 
    #plt.show()

def writeOutput ():
    plt.grid()
    plt.ylim(0,20)
    plt.xlabel(u"x",fontsize=15)
    plt.ylabel(u"y",fontsize=15)
    plt.title(u"Comparação entre as interpolações",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.savefig("output.pdf")

def main():

    points_filename = sys.argv[1]
    interpolate_filename = sys.argv[2]
    interpolate_filename2 = sys.argv[3]
    
    x,y = readPoints(points_filename)
    z1 = readInterpolatePoints(interpolate_filename)
    z2 = readInterpolatePoints(interpolate_filename2)

    showPoints(x,y)
    showInterpolatePoints(z1,"red","Lagrange")
    showInterpolatePoints(z2,"blue","Cplines")

    writeOutput()

if __name__ == "__main__":
    main()

