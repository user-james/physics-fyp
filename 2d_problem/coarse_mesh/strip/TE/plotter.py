from os import system

if __name__ == '__main__':

    command = "gnuplot -e \"set terminal jpeg; set output 'E{:}.jpg'; set  xlabel 'x'; set ylabel 'y'; set dgrid3d 40,40; set hidden3d; splot 'E{:}.txt' u 1:2:3 with lines\"".format(1, 1)
    system(command)
