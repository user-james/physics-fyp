import numpy as np
import matplotlib.pyplot as plt
import sys


def plotter(save, pltname):
    '''
    Uses user input to save figures for given mode(s)
    '''
    while True:
        mode = int(input("Enter mode number you are interested in: "))
        
        if mode < 0:
            break

        data = np.genfromtxt("./data/TE{:d}.txt".format(mode))
        
        x = np.linspace(data[0, 0], data[-1, 0], 500)
        
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax1.plot(data[:, 0], data[:, 1], 'r-', label="Mode")
        ax2.plot(x, 3*np.exp(-x*x*5e10), 'b-', label="Waveguide")
        ax1.set_xlabel("x")
        ax1.set_ylabel("E(x)")
        ax2.set_ylabel("n(x)")
        #plt.plot(-0.000004*np.ones(100), np.linspace(-0.1, 0.1, 100), 'k--')
        #plt.plot(0.000004*np.ones(100), np.linspace(-.1, 0.1, 100), 'k--')
        
        if save:
            ax1.legend(fontsize=11, loc=1)
            ax2.legend(fontsize=11, loc=2)
            plt.savefig("./plots/{:s}TE{:d}_fd.jpg".format(pltname, mode))
            plt.close()
        else:
            ax1.legend(fontsize=11, loc=2)
            ax2.legend(fontsize=11, loc=1)
            plt.show()




if __name__ =='__main__':

    options = ['s', 'p']
    pltname = 'gaussian'

    if (len(sys.argv) > 1) and (sys.argv[1] in options):
        if sys.argv[1] == 's':
            save = True
        else:
            save = False

        plotter(save, pltname)
    else:
        print("No arguments entered")
