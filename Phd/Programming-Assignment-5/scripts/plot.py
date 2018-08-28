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
    plt.ylim(0,10)
    plt.scatter(x,y)

def showInterpolatePoints (z):
    plt.ylim(0,20)
    plt.plot(z[:,0],z[:,1])
    plt.grid() 
    #plt.show()

def writeOutput (x,y,z):
    plt.savefig("output.pdf")

def main():

    points_filename = sys.argv[1]
    interpolate_filename = sys.argv[2]
    
    x,y = readPoints(points_filename)
    z = readInterpolatePoints(interpolate_filename)

    showPoints(x,y)
    showInterpolatePoints(z)

    writeOutput(x,y,z)

if __name__ == "__main__":
    main()

