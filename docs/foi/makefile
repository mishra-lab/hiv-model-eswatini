full:
	make tex && make bib && make tex && make tex

tex:
	pdflatex -interaction=nonstopmode main

bib:
	biber main

clean:
	texclean .

cls:
	cp /home/jesse/dev/tex/pkg/ms/ms/ms.cls .

wc:
	python3.8 ../wc/wc.py abstract intro model.prior model.prop experiment discussion
