set term pdf font "Times New Roman,9"
set output 'time.pdf'
set grid
set multiplot layout 1,1
# set title 'system size: wid*wid, t=1, E=0.2'
# set k left top
# set xtics 0.2
set xlabel 'wid'
set ylabel 'seconds consume for 100 loops'
set xrange [10:130]
plot "res.txt" u 1:2 w lp t "LU","" u 1:3 w lp t "T", '' u 1:5 w lp t 'dosR'
# set yrange [0:10]
unset multiplot
set output
q
