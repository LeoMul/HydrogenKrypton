set t epslatex 9 standalone color size 5.7cm, 3.33cm header '\usepackage{amsmath}' font "Helvetica,10"
set output "titlegraphic.tex"
#set output "cross_sec_book_recreated.pdf"

#set xlabel "Energy / meV"
#set ylabel 'Cross-Sec/ $\rho ^2$'
set xtics format " "
set tics nomirror
set ytics format " "
set border lc "white"
p "cross_sec_Emin0.100Emax3.200N019992h0.001000.dat" u 1:($2*3.57*3.57) w l lw 5 lc rgb "black" t "" 

set output
system("latexmk -pdf titlegraphic.tex")
