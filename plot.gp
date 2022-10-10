set t pdf
set output "cross_sec_book_recreated.pdf"

set xlabel "Energy / meV"
set ylabel "Cross-Sec/ rho^2"

p "cross_sec_e_less_1.dat" w l lc "dark-magenta" t "", "cross_sec_e_more_1.dat" w l lc "dark-magenta" t ""

set output