full:
	make clean && make tex && make bib && make tex && make bib && make tex && make tex

tex:
	pdflatex main

bib:
	bibtex main

clean:
	texclean .

words:
	python3.8 wc/wc.py
