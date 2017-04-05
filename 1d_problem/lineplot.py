import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

if __name__ == '__main__':

    fig, ax = plt.subplots()
    
    ratio_range = np.arange(0.5, 0.0, -.025)
    data1 = np.array([])
    data2 = np.array([])
    data3 = np.array([])
    data4 = np.array([])
    data5 = np.array([])
    for i in ratio_range:
        name = "./TE/{:.3f}_neff.txt".format(i)
        data1 = np.append(data1, np.genfromtxt(name)[0])
        data2 = np.append(data2, np.genfromtxt(name)[1])
        data3 = np.append(data3, np.genfromtxt(name)[2])
        data4 = np.append(data4, np.genfromtxt(name)[3])
        data5 = np.append(data5, np.genfromtxt(name)[4])

    x = np.arange(0.5, 0.0, -.025)*1.6e-5
    

    ax.plot(x, np.sqrt(data1), label=r"$TE_0$")
    ax.plot(x, np.sqrt(data2), label=r"$TE_1$")
    ax.plot(x, np.sqrt(data3), label=r"$TE_2$")
    ax.plot(x, np.sqrt(data4), label=r"$TE_3$")
    ax.plot(x, np.sqrt(data5), label=r"$TE_4$")
    ax.plot(x,3.2*np.ones(len(x)), 'k--')

   
    ax.set_xlabel(r'Waveguide Width')
    ax.set_ylabel(r'$n_{eff}$')
    ax.set_xlim(0.025*1.6e-5, 0.5*1.6e-5)
    ax.set_ylim(3, 3.22)
    plt.legend(loc = 4)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))

    plt.savefig("../important_plots/mode_loss.jpg")
    plt.show()
