export PYTHONPATH := .:$$PYTHONPATH
do.debug:
	python3 model/debug.py
do.model:
	./model/scenario/main.sh
do.scinet:
	sbatch model/scenario/main.sh
do.post.base:
	Rscript post/main.base.r
do.post.foi:
	Rscript post/main.foi.r
do.post.art:
	Rscript post/main.art.r
do.params.distr:
	python3 params/distr.py
do.params.tab:
	python3 params/tab/main.py
do.params.fsw:
	Rscript params/fsw/main.r
do.params.wp:
	Rscript params/wp/main.r
do.profile:
	./profile/run.sh
do.toy:
	python3 toy/Bxph.py
do.tikz:
	cd tikz && make all
