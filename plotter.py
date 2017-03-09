import numpy as np
import matplotlib.pyplot as plt



if __name__=='__main__':

    data = np.genfromtxt("tmvalues1.343061.txt")
    r, a, x, y = data


    c = np.exp(y)*np.cos(x)

    mu, gamma = x/a, y/a

    x1 = np.linspace(-2*a, -a, 100)
    x2 = np.linspace(-a, a, 200)
    x3 = np.linspace(a, 2*a, 100)


    plt.plot(x1, c*np.exp(gamma*x1), 'b-')
    plt.plot(x3, c*np.exp(-gamma*x3), 'b-')
    plt.plot(x2, np.cos(mu*x2), 'r-')

    plt.plot(-a*np.ones(100), np.linspace(-1.5, 1.5, 100), 'k--')
    plt.plot(a*np.ones(100), np.linspace(-1.5, 1.5, 100), 'k--')

    #plt.ylim(0, 1.2)
    plt.savefig("TM1.jpeg")
    plt.show()

