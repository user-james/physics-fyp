import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


def inv_para(x):
    '''
    Function that forms an inverse parabola between [-limit, +limit]
    The function is equal to 1 outside this limit
    '''

    return 3-5e9*x*x

if __name__ == '__main__':

    waveguide = input("Enter name of waveguide function: ")

    data = np.genfromtxt("./effective_indices/{:s}.txt".format(waveguide))

    data_range = np.genfromtxt("./modes/TE0.txt")

    x = np.linspace(data_range[0, 0], data_range[-1, 0], 401)
    xin = x[150:251]
    xleft = x[:150]
    xright = x[251:]

    fig, ax = plt.subplots()
    data = np.sqrt(data)
    count = 0
    for i in data:
        count += 1
        if i < 3 and i > 1:
            if count == 1:
                ax.plot(x, i*np.ones(len(x)), 'r-', label = r"$n_{eff}$")
            else:
                ax.plot(x, i*np.ones(len(x)), 'r-')


    ax.plot(xin, inv_para(xin), 'b-', label="$Waveguide$")
    ax.plot(xleft, np.ones(len(xleft)), 'b-')
    ax.plot(xright, np.ones(len(xright)), 'b-')
    
    save = input("Save file? (y/n) ")

    if save[0] == 'y':
        ax.legend(fontsize = 11)
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$n(x)$')
        ax.set_ylim(0.5, 3.2)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.savefig("./effective_indices/{:s}.jpg".format(waveguide))
        plt.show()
    else:
        plt.legend(fontsize = 11)
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$n(x)$')
        ax.set_ylim(0.5, 3.2)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))
        plt.show()
