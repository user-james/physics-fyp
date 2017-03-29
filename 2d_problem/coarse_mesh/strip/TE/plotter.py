import os
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.pyplot as plt


if __name__=='__main__':

    data = np.genfromtxt("E0.txt")
    x = data[: :40, 0]
    y = data[:40, 1]
    z = data[:, 2]

    x, y = np.meshgrid(x, y)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(x, y, z, rstride=10, cstride=10)

    plt.show()


