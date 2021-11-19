all: opti_book.c
	gcc opti_book.c -O2 -g -lm -o ./exp
	./exp

plot: opti_book.c
	gcc opti_book.c -O3 -lm -o ./exp
	./exp > vals.csv
	gnuplot ploter.p
