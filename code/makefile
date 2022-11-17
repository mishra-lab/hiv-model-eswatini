export PYTHONPATH := .:$$PYTHONPATH
do.debug:
	python3.8 model/debug.py && echo [OK]
do.model:
	python3.8 model/main.py 0 && echo [OK]
do.model.foi:
	python3.8 model/scenario/foi.py && echo [OK]
do.model.art:
	python3.8 model/scenario/art.py && echo [OK]
do.params.distr:
	python3.8 params/distr.py && echo [OK]
do.profile:
	./profile/run.sh && echo [OK]
do.posterior:
	Rscript post/posterior.r && echo [OK]
do.infections:
	Rscript post/infections.r && echo [OK]
do.post.foi:
	Rscript post/foi.r && echo [OK]
do.post.art.1:
	Rscript post/art_obj_1.r && echo [OK]
do.post.art.2:
	Rscript post/art_obj_2.r && echo [OK]
do.tikz:
	cd tikz && make all
do.params.fsw:
	Rscript params/fsw/main.r && echo [OK]
do.params.wp:
	Rscript params/wp/main.r && echo [OK]
do.model.batch:
	python3.8 model/main.py 0 && \
	python3.8 model/main.py 1 && \
	python3.8 model/main.py 2 && \
	python3.8 model/main.py 3 && \
	python3.8 model/main.py 4 && \
	python3.8 model/main.py 5 && \
	python3.8 model/main.py 6 && \
	python3.8 model/main.py 7 && \
	python3.8 model/main.py 8 && \
	python3.8 model/main.py 9 && \
	echo [OK]
