# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt

def read_solution (filename):
    data = np.genfromtxt(filename, delimiter=' ')
    return data

def show_solution (data):
    x = data[:,0]
    y_aprox = data[:,1]
    y_analit = data[:,2]

    plt.plot(x,y_aprox,label="aprox")
    plt.plot(x,y_analit,label="analit")

def write_solution ():
    #plt.grid()
    plt.xlabel(u"x",fontsize=15)
    plt.ylabel(u"u",fontsize=15)
    plt.title(u"Analitical x Aproximation",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.show()
    
def main():

    if (len(sys.argv) != 2):
        print("==================================================================")
        print("Usage:> python plot.py <solution_file>")
        print("==================================================================")
        print("<solution_file> = File with the solution from the simulation")
        print("==================================================================")
        sys.exit(1)

    solution_filename = sys.argv[1]
    
    data = read_solution(solution_filename)

    show_solution(data)
    write_solution()


if __name__ == "__main__":
    main()
