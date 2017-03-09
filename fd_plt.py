import numpy as np
import matplotlib.pyplot as plt
import sys


def plotter(save):
    '''
    Uses user input to save figures for given mode(s)
    '''
    while True:
        mode = int(input("Enter mode number you are interested in: "))
        
        if mode < 0:
            break

        data = np.genfromtxt("./data/TE{:d}.txt".format(mode))


        plt.plot(data[:, 0], data[:, 1])
        plt.plot(-0.000004*np.ones(100), np.linspace(-0.1, 0.1, 100), 'k--')
        plt.plot(0.000004*np.ones(100), np.linspace(-.1, 0.1, 100), 'k--')
        
        if save:
            plt.savefig("./plots/TE{:d}_fd.jpg".format(mode))
            plt.close()
        else:
            plt.show()




if __name__ =='__main__':

    options = ['s', 'p']
    
    if (len(sys.argv) > 1) and (sys.argv[1] in options):
        if sys.argv[1] == 's':
            save = True
        else:
            save = False

    plotter(save)
