do.all:
	make do.fig && make do.abstract
do.abstract:
	pdflatex -shell-escape main
do.fig:
	cd fig && pdflatex main && convert -density 600 -flatten main.pdf main.png