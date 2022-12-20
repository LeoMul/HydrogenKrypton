set t epslatex 9 standalone color size 8.6cm, 5cm header '\usepackage{amsmath}' font "Helvetica,10"
set output "plot.tex"
#set output "cross_sec_book_recreated.pdf"

set xlabel '$E$/ meV'
set ylabel '$\sigma$/ Å$^2$'

p "cross_sec_Emin0.250Emax3.500N004992.dat" u 1:($2*3.57*3.57)w l lc "black" t ""

set output
system("latexmk -pdf plot.tex")
