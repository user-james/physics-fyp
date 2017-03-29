#!/usr/bin/gnuplot

set dgrid3d 40, 40
set hidden3d
splot "E0.txt" u 1:2:3 with lines
pause -1
