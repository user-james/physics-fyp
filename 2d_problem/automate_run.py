from subprocess import call
import numpy as np


if __name__ == '__main__':

    scales_init = np.arange(0.6, 0.4, -.05)
    scales_end = np.arange(0.4, 0.05, -0.01) 

    scales = np.hstack((scales_init, scales_end))
    for scale in scales:
        command = "{:.2f}".format(scale)
        call(["./f_diff_2d", "s", command])
