import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':

    x = np.linspace(0, 1, 100)

    plt.plot(x, x, 'r-', label=r'$n_{eff}$')
    plt.plot(x, 8*x, 'b-', label=r'$Waveguide$')

    plt.legend(fontsize = 20)
    plt.show()

