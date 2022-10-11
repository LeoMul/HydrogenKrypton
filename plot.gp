set t epslatex 9 standalone color size 8.6cm, 5cm header '\usepackage{amsmath}' font "Helvetica,10"
set output "plot.tex"
#set output "cross_sec_book_recreated.pdf"

set xlabel "Energy / meV"
set ylabel 'Cross-Sec/ $\rho ^2$'

p "cross_sec_e_less_1.dat" w l lc "dark-magenta" t "", "cross_sec_e_more_1.dat" w l lc "dark-magenta" t ""

set output
system("latexmk -pdf plot.tex")
