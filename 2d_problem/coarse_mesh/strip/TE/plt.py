import os 


os.system("gnuplot -p")
os.system("set dgrid3d 40, 40")
os.system("set hidden3d")
os.system("splot 'E0.txt' u 1:2:3 with lines")

