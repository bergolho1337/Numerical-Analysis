import numpy as np
import matplotlib.pyplot as plt

def main():

    #x = np.linspace(-1, 1, 100)
    #y = 0.5*x + 0.2*np.random.rand(100)
    x = np.array([0.0, 1.0, 2.0, 6.0, 7.0])
    y = np.array([0.0, 0.0, 1.0, 2.0, 3.0])

    #np.polynomial.chebyshev.Chebyshev.fit
    p = np.polynomial.Chebyshev.fit(x,y,1)
    plt.scatter(x,y)
    plt.plot(x,p(x),c="red")
    plt.show()


if __name__ == "__main__":
    main()