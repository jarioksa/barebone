all: basicmethods.pdf
barebone.tex: barebone.Rnw
	"$(R_HOME)/bin/R" CMD Sweave barebone.Rnw
basicmethods.pdf: basicmethods.tex barebone.tex
	"$(R_HOME)/bin/R" CMD texi2pdf --clean basicmethods.tex
	-rm -f barebone-0*.* Rplots.*
