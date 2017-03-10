import numpy as np
import matplotlib.pyplot as plt


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

    count = 0
    for i in data:
        count += 1
        if i < 3:
            if count == len(data):
                plt.plot(x, i*np.ones(len(x)), 'r-', label = r"$n_{eff}$")
            else:
                plt.plot(x, i*np.ones(len(x)), 'r-')


    plt.plot(xin, inv_para(xin), 'b-', label="$Waveguide$")
    plt.plot(xleft, np.ones(len(xleft)), 'b-')
    plt.plot(xright, np.ones(len(xright)), 'b-')
    
    save = input("Save file? (y/n) ")

    if save[0] == 'y':
        plt.legend(fontsize = 11)
        plt.savefig("./effective_indices/{:s}.jpg".format(waveguide))
        plt.show()
    else:
        plt.legend(fontsize = 11)
        plt.show()
