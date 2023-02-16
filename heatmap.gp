set t epslatex 9 standalone color size 5in,3in header '\usepackage{amsmath}' font "Helvetica,12"
set output "heatmap_analy.tex"

#set pm3d interpolate 2,2


set view map
#set dgrid3d
#set xtics 0,1,4
set xlabel '$E$ / meV'
set ytics ("0" 0, '$\pi/4$' 0.785, '$\pi/2$' 1.57, '$3\pi/4$' 2.356, '$\pi$' 3.14)
set yrange[0:3.1459]
set xrange[0.3:3.2]
set ylabel '$\theta$' rotate by 0  
set cblabel '\Large $\frac{d \sigma}{ d \Omega}$' rotate by 0 
#set size 1,1
#set cbrange [0:0.37802209909052564]
maxcb = 1000
set cbrange[0:maxcb]
set palette defined ( 0 "#000000", maxcb/12*0.4 "#7e386c", maxcb/12 * 2*0.5 "#a3597c", maxcb/12*3*0.6 "#af5a6c", maxcb/12*4*0.7 "#b88260", maxcb/12*5*0.8 "#b49360", maxcb/12*6*0.9 "#a9b370", maxcb/12*7*0.95 "#9ec07e", maxcb/12*8 "#a3d6ae", maxcb/12*9 "#addedb", maxcb/12*10 "#ccdef6", maxcb/12*11 "#f7e9fa", maxcb/12*12 "#ffffff")


set rmargin screen 0.8
set lmargin screen 0.15
#set palette rgbformulae 22,13,10
splot "dcsa_data.dat" u 2:1:3 t "" with image

set output
system("latexmk -pdf heatmap_analy.tex")
