all: basicmethods.pdf
basicmethods.pdf: basicmethods.tex barebone.tex
	texi2dvi --pdf --clean basicmethods.tex
	-rm -f barebone-0*.* Rplots.*
