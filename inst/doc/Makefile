all: basicmethods.pdf
barebone.tex: barebone.Rnw
	"$(R_HOME)/bin/R" CMD Sweave barebone.Rnw
basicmethods.pdf: basicmethods.tex barebone.tex
	texi2dvi --pdf --clean basicmethods.tex
	-rm -f barebone-0*.* Rplots.*
