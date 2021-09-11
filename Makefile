CFLAGS=-W -Wall -O3 -ffast-math -march=native -funroll-loops -fmodulo-sched -finline-limit=100000

all: fec

OBJECT= fec.o derivs.o ode-rk-4.o ode-rkf-45.o get-param.o diagonalize-sym-3.o reset.o

fec: ${OBJECT}
	${CC} ${CFLAGS} ${OBJECT} -o fec -lm

USAGE.pdf: USAGE.tex
	latex USAGE
	latex USAGE
	dvips USAGE -o
	ps2pdf USAGE.ps

fec.exe:
	make clean
	sed -e s/-march=native// Makefile | make CC=mingw32-gcc -f -
	mv fec fec.exe
	make clean

tar:
	make veryclean
	make USAGE.pdf
	make fec.exe
	BASE=`basename $$PWD` && cd .. && tar cvfz $${BASE}.tar.gz $${BASE}
	BASE=`basename $$PWD` && cd .. && rm -f $${BASE}.zip && zip -r $${BASE}.zip $${BASE}

clean:
	rm -f fec *.o *~ *.aux *.dvi *.log *core s.out s-out.eps *.bak *.png *.ps

veryclean: clean
	rm -f *.pdf *.exe
