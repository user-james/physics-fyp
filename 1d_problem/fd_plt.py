import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.ticker import FormatStrFormatter

def plotter(save, pltname):
    '''
    Uses user input to save figures for given mode(s)
    '''
    width = 1.6e-5
    while True:
        mode = int(input("Enter mode number you are interested in: "))
        if mode < 0:
            break
        
        ratio = float(input("Ratio = "))
        right = width*ratio/2
        left = -right

        data = np.genfromtxt("./TE/{:.3f}_E{:d}.txt".format(ratio, mode))
        
        x = np.linspace(data[0,0], data[-1,0], 801)
        xin = x[300:501]
        xleft = x[:300]
        xright = x[501:]


        fig, ax1 = plt.subplots()
        #ax2 = ax1.twinx()
        ax1.plot(data[:, 0], data[:, 1], 'r-', label="Mode")
        up, down = max(data[:, 1]), min(data[:, 1])
        yspace = np.linspace(down, up, 100)
        #ax2.plot(xin, inv_para(xin), 'b-', label="Waveguide")
        #ax2.plot(xleft, np.ones(300), 'b-')
        #ax2.plot(xright, np.ones(300), 'b-')
        ax1.set_xlabel(r"$x$")
        ax1.set_ylabel(r"$E(x)$")
        #ax1.set_xlim(-4e-6, 4e-6)
        #ax2.set_ylabel("n(x)")
        ax1.plot(left*np.ones(100), yspace, 'k--')
        ax1.plot(right*np.ones(100), yspace, 'k--')
        
        if save:
            #ax2.legend(fontsize=11, loc=1)
            ax1.legend(fontsize=11, loc=2)
            ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0e'))
            plt.savefig("../important_plots/TE{:d}_fd.jpg".format(mode))
            plt.close()
        else:
            ax1.legend(fontsize=11, loc=2)
            ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0e'))
            #ax2.legend(fontsize=11, loc=1)
            plt.show()

def inv_para(x):
    '''
    Function that forms an inverse parabola between [-limit, +limit]
    The function is equal to 1 outside this limit
    '''

    return 3-5e9*x*x


if __name__ =='__main__':

    options = ['s', 'p']
    pltname = 'inv_parabola'

    if (len(sys.argv) > 1) and (sys.argv[1] in options):
        if sys.argv[1] == 's':
            save = True
        else:
            save = False

        plotter(save, pltname)
    else:
        print("No arguments entered")
