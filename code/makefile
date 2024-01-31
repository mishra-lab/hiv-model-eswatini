export PYTHONPATH := .:$$PYTHONPATH
do.debug:
	python3.8 model/debug.py && echo [OK]
do.model.imis:
	python3.8 model/scenario/imis.py case=base b=0 && echo [OK]
do.model.foi:
	python3.8 model/scenario/foi.py && echo [OK]
do.model.art:
	python3.8 model/scenario/art.py && \
	# python3.8 model/scenario/art.py b=0 && \
	# python3.8 model/scenario/art.py b=1 && \
	# python3.8 model/scenario/art.py b=2 && \
	# python3.8 model/scenario/art.py b=3 && \
	echo [OK]
do.post.base:
	Rscript post/main.base.r && echo [OK]
do.post.foi:
	Rscript post/main.foi.r && echo [OK]
do.post.art:
	Rscript post/main.art.r && echo [OK]
do.params.distr:
	python3.8 params/distr.py && echo [OK]
do.params.tab:
	python3.8 params/tab/main.py && echo [OK]
do.params.fsw:
	Rscript params/fsw/main.r && echo [OK]
do.params.wp:
	Rscript params/wp/main.r && echo [OK]
do.profile:
	./profile/run.sh && echo [OK]
do.toy:
	python3.8 toy/Bxph.py && echo [OK]
do.tikz:
	cd tikz && make all
