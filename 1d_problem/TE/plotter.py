import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':

    data = np.genfromtxt("0.200_E0.txt")
    
    fig, ax = plt.subplots()

    ax.plot(data[:, 0], data[:, 1])
    
    ax.set_ylim(0, 0.12)
    ax.set_xlim(-4e-6, 4e-6)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$E_y$')
    plt.show()
