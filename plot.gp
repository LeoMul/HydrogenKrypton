set t epslatex 9 standalone color size 8.6cm, 5cm header '\usepackage{amsmath}' font "Helvetica,10"
set output "plot.tex"
#set output "cross_sec_book_recreated.pdf"

set xlabel '$E$/ meV'
set ylabel '$\sigma$/ Å$^2$'

p "cross_sec_Emin0.100Emax0.300N000048.dat" u 1:($2*3.57*3.57) w l lc rgb  "black" t "" , "cross_sec_Emin0.300Emax1.000N000996.dat" u 1:($2*3.57*3.57)w l lc "black" t "","cross_sec_Emin1.000Emax3.200N001500.dat" u 1:($2*3.57*3.57)w l lc "black" t ""
set output
system("latexmk -pdf plot.tex")
