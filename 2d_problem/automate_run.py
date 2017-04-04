from subprocess import call
import numpy as np


if __name__ == '__main__':

    scales_init = np.arange(1.0, 0.6, -.05)
    #scales_end = np.arange(0.4, 0.05, -0.01) 

    #scales = np.hstack((scales_init, scales_end))
    total = len(scales_init)
    count = 0
    for scale in scales_init:
        count += 1
        command = "{:.2f}".format(scale)
        call(["./f_diff_2d", "s", command])
        print("\n{:d} iteration done our of {:d}\n".format(count, total))
