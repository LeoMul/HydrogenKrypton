set t epslatex 9 standalone color size 6.88cm, 4cm header '\usepackage{amsmath}' font "Helvetica,12"
set output "lj.tex"

set xrange [0.5:2.5]
set yrange [-1.3:5.1]
set xlabel ('$r/\rho$')
set ylabel ('$V_{\text{LJ}}/\epsilon$')
set ytics -1,2,5
f(x) = (1/x)**12 - 2*((1/x)**6)
set samples 1000
p f(x) w l lc "black" lw 5 t ""

set output

system("latexmk -pdf lj.tex")
