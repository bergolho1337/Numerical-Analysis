import sys
import numpy as np
import matplotlib.pyplot as plt

def readPoints (filename):
    filename = sys.argv[1]
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
    plt.scatter(x,y)
    plt.ylim(0,30)
    plt.show()



def main():

    filename = sys.argv[1]
    
    x,y = readPoints(filename)
    showPoints(x,y)

if __name__ == "__main__":
    main()

