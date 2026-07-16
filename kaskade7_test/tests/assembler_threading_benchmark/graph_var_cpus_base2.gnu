set style line  7 lt 1 lc rgb "red" pt 0 lw 4
set style line  8 lt 2 lc rgb "red" pt 0 lw 3
set style line  9 lt 3 lc rgb "red" pt 0 lw 3
set style line 10 lt 4 lc rgb "red" pt 0 lw 3
set style line  6 lt 5 lc rgb "red" pt 0
set style line  5 lt 6 lc rgb "red" pt 0
set style line  4 lt 7 lc rgb "red" pt 0

set style line 11 lt 1 lc rgb "blue" pt 0 lw 4
set style line 12 lt 2 lc rgb "blue" pt 0 lw 3
set style line 13 lt 3 lc rgb "blue" pt 0 lw 3
set style line 14 lt 4 lc rgb "blue" pt 0 lw 3
set style line 15 lt 5 lc rgb "blue" pt 0 lw 3
set style line 16 lt 6 lc rgb "blue" pt 0

set style line 17 lt 1 lc rgb "sea-green" pt 0 lw 4
set style line 18 lt 2 lc rgb "sea-green" pt 0 lw 3

set style line 19 lt 1 lc rgb "black" pt 0 lw 4
set style line 20 lt 2 lc rgb "black" pt 0 lw 3
set style line 21 lt 3 lc rgb "black" pt 0 lw 3

set style line 24 lt 1 lc rgb "yellow" pt 0 lw 4
set style line 25 lt 2 lc rgb "yellow" pt 0 lw 3


#  plots for multithreading in Kaskade7.2 assembly: stiffness matrix (2D)
getenv(var) = system(sprintf("echo $%s", var))
maxthreads = getenv("maxthreads")
set output "speedup_vs_threads_2d_first".maxthreads."_base2.eps"
set terminal postscript enhanced eps 22 color
set title "2D assembly (stiffness matrix) ( graph: 2*threads(2)/threads(n) ) \n on ".getenv("HOSTNAME")
set xrange[1:@maxthreads]
set yrange[0:@maxthreads]
set xlabel "#threads"
set ylabel "speedup"
set key top left spacing 1
set pointsize 0.5
plot x  t 'optimum' with linespoints ls 8, \
"assembsummary_order6.stat" using 1:3 t 'p=6' with linespoints ls 11,\
"assembsummary_order5.stat" using 1:3 t 'p=5' with linespoints ls 7,\
"assembsummary_order4.stat" using 1:3 t 'p=4' with linespoints ls 19,\
"assembsummary_order3.stat" using 1:3 t 'p=3' with linespoints ls 24,\
"assembsummary_order2.stat" using 1:3 t 'p=2' with linespoints ls 17,\
"assembsummary_order1.stat" using 1:3 t 'p=1' with linespoints ls 1

