set t epslatex 9 standalone color size 8.6cm, 5cm header '\usepackage{amsmath}' font "Helvetica,12"
set output "plot_with_labels_and_exp.tex"
#set output "cross_sec_book_recreated.pdf"

set xlabel '$E$/ meV'
set ylabel '$\sigma$/ Å$^2$'
set ytics 0,140,700
set label '$l = 4$' at 0.6,600
set label '$l = 5$' at 1.5,450
set label '$l = 6$' at 2.5,405
set key samplen 2

scale = 8.5

set xrange[0:3.2]
set yrange[0:700]

p "cross_sec_Emin0.100Emax0.300N000048.dat" u 1:($2*3.57*3.57) w l lw 5 lc rgb "black" t '\small Numerical' , "cross_sec_Emin0.300Emax1.000N000996.dat" u 1:($2*3.57*3.57)w l lw 5 lc "black" t  "","cross_sec_Emin1.000Emax3.200N001500.dat" u 1:($2*3.57*3.57)w l lc "black" lw 5 t "", "toenniesexperimentaldata.txt" u 1:($2*scale) ps 1.5 pt 2 lc rgb "red" t '\small Experiment'
set output
system("latexmk -pdf plot_with_labels_and_exp.tex")
