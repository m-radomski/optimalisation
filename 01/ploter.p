set terminal png medium size 1280,960
set output "plot.png"

by3(x) = (((int(x)%3)+1)/6.)
by4(x) = (((int(x)%4)+1)/7.)
rgbf(x) = x*51*32768 + (11-x)*51*128 + int(abs(5.5-x)*510/9.)
rgb(r,g,b) = 65536 * int(rgbf(r)) + 256 * int(rgbf(g)) + int(rgbf(b))

set datafile separator ","
set samples 1000, 1000
set xrange [-150:150]
set yrange [-2:2]

f(x)=-cos(0.1 * x)*exp(-((0.1*x - 2*pi)**2)) + 0.002 * (0.1*x)**2
plot f(x), "vals.csv" u ($3):(abs($4/100)):($4-$3):(0.0):(rgb($1,$2,$3)) w xerrorbars lc rgb variable, \
"vals.csv" u (x=$5):(y=f($5)) w p pt 2 ps 3 lw 3 lc rgb "red" title "Found minima [Fibonacci]", \
"vals.csv" u (x=$6):(y=f($6)) w p pt 7 ps 2 lc rgb "green" title "Found minima [Lagrange]"
