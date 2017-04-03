import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter



if __name__=='__main__':

    data = np.genfromtxt("./modes/values1.355118.txt")
    r, a, x, y = data


    c = np.exp(y)*np.cos(x)
    print(c)
    mu, gamma = x/a, y/a

    x1 = np.linspace(-2*a, -a, 100)
    x2 = np.linspace(-a, a, 200)
    x3 = np.linspace(a, 2*a, 100)
    
    fig, ax = plt.subplots()

    ax.plot(x1, c*np.exp(gamma*x1), 'b-')
    ax.plot(x3, c*np.exp(-gamma*x3), 'b-')
    ax.plot(x2, np.cos(mu*x2), 'r-')

    ax.plot(-a*np.ones(100), np.linspace(-1.5, 1.5, 100), 'k--')
    ax.plot(a*np.ones(100), np.linspace(-1.5, 1.5, 100), 'k--')

    ax.set_ylim(0, 1.2)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$E_y$')
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    plt.savefig("./plots/TE0.jpeg")
    plt.show()

