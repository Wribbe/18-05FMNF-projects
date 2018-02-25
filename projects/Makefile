FIGS := figs/*.pdf
FIGS += figs/*.png

all: p1.pdf p2.pdf

nopy:
	pdflatex p2.tex

p2.pdf : p2.tex p2.py ${FIGS}
	python p2.py
	pdflatex $^

p1.pdf : p1.tex ${FIGS}
	pdflatex $^
