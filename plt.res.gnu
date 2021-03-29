set term pdf font "Times New Roman,9" ps 0.6
set output 'time_1subband.pdf'
set grid
set multiplot layout 1,1
# set title 'system size: wid*wid, t=1, E=0.2'
# set k left top
# set xtics 0.2
set xlabel 'wid'
set ylabel 'seconds consume for 100 loops'
set xrange [40:140]
set yrange [0:7.5]
plot "res.1subband.txt" u 1:2 w lp t "LU","" u 1:3 w lp t "T", '' u 1:5 w lp t 'dosR'
unset multiplot
set output
q
