set terminal postscript eps enhanced "Times" 18

set grid  xtics lt 4 lc rgbcolor "grey"
set grid  ytics lt 4 lc rgbcolor "grey"
set grid mxtics lt 4 lc rgbcolor "grey"
set grid mytics lt 4 lc rgbcolor "grey"

set size ratio 0.5

set output "fig_SSP_crweno5_4_MaxDtStages.eps"
set xlabel "Number of stages"
set ylabel "Stage CFL"
set xrange [3:14]
set yrange [0.4:0.7]
set nokey
plot 'MaxDtStages_4thorder.dat' u 1:2 w lp pt 7 ps 2 lt 4 lw 2 lc rgbcolor "blue"
