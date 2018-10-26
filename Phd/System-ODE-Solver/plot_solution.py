import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def plot_solution (data):
    plt.grid()
    plt.plot(data[:,0],data[:,1],label="aprox",c="blue",linestyle='--')
    plt.plot(data[:,0],data[:,3],label="analit",c="red")
    plt.xlabel("t",fontsize=15)
    plt.ylabel("y",fontsize=15)
    plt.title("Analitical x Aproximation",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.show()    

def main():

    if (len(sys.argv) != 2):
        print("==========================================================")
        print("Usage:> %s <input_file>" % (sys.argv[0]))
        print("==========================================================")
        sys.exit(1)
    else:
        input_file = sys.argv[1]

        data = np.genfromtxt(input_file)

        plot_solution(data)        
    

if __name__ == "__main__":
    main()